#' tidysurv package
#'
#' @import lazyeval
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import broom
#' @import stringr
#' @import flexsurv
#' @import survival
#'
#' @importFrom tibble rownames_to_column
#' @importFrom broom tidy
#' @importFrom purrr map map2 pmap map_dbl map_chr map_lgl
#'             map2_chr map2_dbl map2_lgl transpose flatten_chr flatten_dbl flatten_lgl
#'             walk walk2
#' @importFrom graphics plot
#' @importFrom stats sd terms update integrate delete.response
#'             as.formula coef fitted glm lm median na.omit poly
#'             predict var quantile model.response model.frame
#'             na.pass na.exclude na.fail model.matrix model.weights
#'             .getXlevels .checkMFClasses reformulate
#' @importFrom grDevices dev.off pdf
#'
#' @docType package
#' @name tidysurv
NULL
#> NULL

#' Make a tidy survival-object
#'
#' See \code{?tidysurv.formula} and \code{tidysurv.cr_survreg} for more details.
#'
#' @param object Current methods support a formula or an object of type \code{cr_survreg}
#' @param ... Further arguments to be passed to methods.
#'
#' @return A tidy dataframe, describing survival over time
#' @export
#'
tidysurv <- function(object, ...) {
  UseMethod('tidysurv')
}

#' Make a tidy survfit object from a formula
#'
#' @param object A formula, as would be specified for \code{survival::survfit.formula}
#' @param data A dataframe
#' @param na.action a missing-data filter function
#' @param time_period You can choose to bin events into time-periods. This can be helpful if you want
#'        plot the rate within a time-bin, rather than the cumulative rate of incidence.
#' @param max_num_levels Since this function creates a separate survival-curve for each distinct
#'   level of predictors (the right hand side of the \code{formula}), it can sometimes be the case
#'   that you'll accidentally create hundreds or thousands of distinct survival curves, because you
#'   accidentally passed the wrong column-name to the formula. This argument ensures that you don't
#'   accidentally eat up time computing these curves when you didnt mean to. Default is 200.
#' @param ... Further arguments to be passed to \code{survival::survfit.formula}
#'
#' @return A tidy dataframe
#' @export
#'
tidysurv.formula <- function(object, data, na.action = na.exclude, time_period = NULL, max_num_levels = 200, ...) {
  stopifnot( is.data.frame(data) )
  if( is.element('..strata',colnames(data)) )
    stop(call. = FALSE, "The column name '..strata' is reserved, please rename.")

  if (!is.null(time_period)) {
    time_col_name <- as.character(lazyeval::as.lazy(lazyeval::f_lhs(object), env = parent.frame())$expr$time)
    if (!is.element(time_col_name, colnames(data))) {
      stop(call. = FALSE,
           "For the `time_period` arg to work, your `Surv` object in your formula should reference a column in your dataset. Instead, found: ", time_col_name)
    }
    time_grid <- seq(from = time_period, to = max(data[[time_col_name]],na.rm=TRUE), by = time_period)
    time_mat <- matrix(data = rep(time_grid, each = nrow(data)), nrow = nrow(data))
    data[[time_col_name]] <- matrixStats::rowCollapse(x = time_mat, idxs = apply(X = abs(data[[time_col_name]] - time_mat), FUN = which.min, MARGIN = 1))
  }

  ## collect covariates:
  tobj <- stats::delete.response(stats::terms(object))
  if ( length(all.vars(tobj)) > 0 ) {
    mf <- stats::model.frame(formula = tobj, data=data, na.action=na.pass)
    mfd <- dplyr::distinct(mf)
    strata_names <- colnames(mfd)
    mfd$..strata <- as.character(1:nrow(mfd))
    if (nrow(mfd) > max_num_levels)
      stop(call. = FALSE, "The right-hand-side of the formula defines more than ",
           max_num_levels,
           " unique levels. Increase `max_num_levels` if this was intentional.")
    newdata <- bind_cols( data[,all.vars(update(object, .~1)),drop=FALSE],
                          dplyr::left_join(x = mf, y = mfd, by = colnames(mf)) )
    sfit <- survival::survfit(formula = update(object, .~..strata), data = newdata, na.action = na.action, ... = ...)
  } else {
    sfit <- survival::survfit(formula = object, data = data, na.action = na.action, ... = ...)
    strata_names <- NULL
  }

  ## convert strata -> covariates:
  df_tidy <- tidy.survfit(sfit)
  if (inherits(sfit, "survfitms"))
    df_tidy <- dplyr::rename(.data = df_tidy, .surv_state = state)
  if ( length(all.vars(tobj)) > 0 )  {
    df_tidy$..strata <- stringr::str_match(string = df_tidy$strata, pattern = "\\.\\.strata=(.*)")[,2]
    out <- dplyr::left_join(x = df_tidy, y = mfd, by = '..strata')
    out <- dplyr::as_data_frame(dplyr::select(.data = out, -`..strata`, -strata))
  } else {
    out <- dplyr::as_data_frame(df_tidy)
  }

  if (!is.null(time_period)) {
    get_binom_ci <- failwith(f = function(x,n) broom::tidy(binom.test(x = sum(x), n = sum(n) )),
                             default = data_frame(estimate=NA), quiet = F)

    df_rate_ci <- out %>%
      rowwise() %>%
      do( get_binom_ci(x = round(.$n.event), n = round(.$n.risk) )) %>%
      select(rate = estimate, rate.high = conf.high, rate.low = conf.low)
    out <- bind_cols(out, df_rate_ci)
  }

  ## out
  class(out) <- c('tidysurv',class(out))
  attr(out, 'tidysurv') <- list(strata_names = strata_names, method = 'formula', time_period = time_period)
  return( out )

}

#' Plot a tidysurv object
#'
#' @param x A tidysurv object
#' @param mapping An aesthetic mapping for ggplot, created by one of the \code{aes} functions.
#' @param states For a competing risks model: which states should be included in the graph?
#' @param type What type of plot? Cumulative/survival, or rate? If the latter, need to have set
#'             \code{time_period} in \code{tidysurv}.
#' @param ... Further arguments to be passed to ggplot
#'
#' @return A ggplot object, with lines for each survival curve, and ribbons for confidence bands.
#'   Defaults to a separate survival curve (with a distinct color) for each unique combination of
#'   predictors that was passed to \code{tidysurv}. Also grouped by `.surv_state` for competing risk models,
#'   with the linetype indicating which state.
#' @export
plot.tidysurv <- function(x, mapping = NULL, type = 'cumulative', states = NULL, ...) {

  if ('.surv_state' %in% colnames(x)) {
    possible_states <- x$.surv_state
    if (is.null(states))
      states <- setdiff(possible_states, '')
    x <- dplyr::filter(x, .surv_state %in% states)
    class(x) <- c('tidysurv', class(x))
  }

  type <- match.arg(arg = type, choices = c('cumulative','rate'))

  if (type == 'cumulative') {
    the_ribbon <- geom_ribbon(aes(ymin=conf.low, ymax=conf.high), color=NA, alpha=.10)
  } else {
    the_ribbon <- geom_ribbon(aes(ymin=rate.low, ymax=rate.high), color=NA, alpha=.10)
  }

  g <- ggplot(x, mapping = mapping, type = type, ... = ...) +
    geom_line() +
    the_ribbon
  if (attr(x,'tidysurv')$method == 'cr_survreg') {
    if (type == 'rate')
      g <- g + geom_line(mapping = aes(y = fitted_hazard), size=1.2,alpha=.75)
    else
      g <- g + geom_line(mapping = aes(y = fitted), size=1.2,alpha=.75)
  }

  return(g)
}

#' Ggplot method for tidysurv objects
#'
#' @param data A tidysurv object
#' @param mapping An aesthetic mapping for ggplot, created by one of the \code{aes} functions.
#' @param type What type of plot? Cumulative/survival, or rate? If the latter, need to have set
#'             \code{time_period} in \code{tidysurv}.
#' @param ... Further arguments to be passed to ggplot
#' @param environment See \code{?ggplot2::ggplot}
#'
#' @return A ggplot object without any layers, but with sensible defaults for aesthetic mappings:
#'   each distinct predictor in a group (with distinct colors), and (if a competing risk model) each
#'   state in a separate group (with distinct line-types).
#' @export
ggplot.tidysurv <- function(data, mapping = NULL, type = 'cumulative', ..., environment = parent.frame()) {

  strata_names <- attr(data,'tidysurv')$strata_names

  if ( length(strata_names) > 1 ) {
    data$..group <- do.call(interaction, data[ , strata_names, drop=FALSE])
    group_name <- '..group'
  } else if (length(strata_names)==1) {
    group_name <- paste0("`", strata_names, "`")
  } else {
    data$`1` <- 1
    group_name <- '1'
  }

  type <- match.arg(arg = type, choices = c('cumulative','rate'))

  if (type != 'rate') {
    the_ribbon <- geom_ribbon(aes(ymin=conf.low, ymax=conf.high), color=NA, alpha=.10)
    y_aes <- ~estimate
  } else {
    if (is.null(attr(data,'tidysurv')$time_period))
      stop(call. = FALSE, "If `type= 'rate'`, then a this tidysurv object needs to have a `time_period.` See `?tidysurv`.")
    y_aes <- ~rate
  }

  if ('.surv_state' %in% colnames(data) && dplyr::n_distinct(data[['.surv_state']]) > 1) {
    default_aes <- aes_(x = ~time, y =y_aes,
                        group = as.formula(paste0('~interaction(.surv_state,',group_name, ')')),
                        color = reformulate(group_name),
                        linetype = ~.surv_state)
  } else {
    default_aes <- aes_(x = ~time, y = y_aes,
                        group = reformulate(group_name),
                        color = reformulate(group_name))
  }

  mapping <- update_aes(default_aes, mapping)
  class(data) <- setdiff(class(data),'tidysurv')
  ggplot(data = data, mapping = mapping, ... = ..., environment = environment)
}

#' Update ggplot aes
#'
#' @param old_aes The old aes
#' @param new_aes New aes to merge in. Replaces any conflicts in old_aes
#'
#' @return An updated aes for ggplot mapping argument.
#'
update_aes <- function(old_aes, new_aes) {
  if (is.null(new_aes)) {
    return(old_aes)
  } else {
    old_aes <- old_aes[setdiff(names(old_aes), names(new_aes))]
    mapping_arg <- structure(c(old_aes, new_aes), class='uneval')
    return(mapping_arg)
  }
}


#' Easy interface for parametric competing-risks survival-regression
#'
#' @param time_col_name A character indicating the column-name of the 'time' column.
#' @param event_col_name A character indicating the column-name for events. The values in this
#'   column (either factor or character), should match the names of the \code{list_of_list_of_args},
#'   except for the value indicating censoring.
#' @param data A dataframe.
#' @param list_of_list_of_args A list, whose names correspond to each possible type of event. These
#'   should correspond to the values in the \code{event_col_name} column. The elements of this list
#'   are themselves lists— arguments to be passed to the survival-regression modelling function (as
#'   specified by \code{method}— currently only \code{flexsurvreg} is supported).
#' @param time_lb_col_name Optional. A character indicating the column-name of the time
#'   'lower-bound'. See the vignette for details.
#' @param time_start_col_name Optional. A character indicating the column-name of the start
#'   (truncation) times. See the vignette for details.
#' @param method A character string naming the function to be called for survival-regression
#'   modelling. Currently only supports \code{flexsurvreg}.
#' @param memoize When this function runs, it calls the fitting method (e.g., \code{flexsurvreg}),
#' once per event-type. If you're only tweaking the model for one event-type, you might not want to
#' waste time re-fitting the model for all of the other event types. When \code{memoize} is TRUE (the
#' default), calls to the fitting method are cached: if called with the exact same arguments, the
#' code isn't re-run. This means that you can tweak the model for one event-type at a time, without
#' having to worry about re-fitting the other models on each function call.
#'
#' @return An object of type \code{cr_survreg}, with plot and summary methods.
#' @export
#'
cr_survreg <- function(time_col_name,
                       event_col_name,
                       data,
                       list_of_list_of_args,
                       time_lb_col_name = NULL,
                       time_start_col_name = NULL,
                       method = 'flexsurvreg',
                       memoize = TRUE) {


  ## memoized functions:
  mem_flexsurvreg <- memoise::memoise(function(formula, anc=NULL, data, weights=NULL,bhazard=NULL,subset=NULL,
                                               na.action= na.exclude, dist, inits=NULL, fixedpars=NULL, dfns=NULL,
                                               aux=NULL, cl=.95, integ.opts = NULL, sr.control = survival::survreg.control(), ...) {
    if ( (!is.null(weights)) | (!is.null(bhazard)) | (!is.null(subset)) | (!is.null(integ.opts)) )
      stop("Weights, bhazard, integ.opts, and subset are non supported if `memoize= TRUE`.")
    if (is.null(inits))
      flexsurv::flexsurvreg(formula = formula, anc=anc, data = data,
                            na.action=na.action, dist=dist, fixedpars=fixedpars,
                            dfns=dfns, aux=aux, cl=cl, sr.control = sr.control, ... = ...)
    else
      flexsurv::flexsurvreg(formula = formula, anc=anc, data = data, inits=inits,
                            na.action=na.action, dist=dist, fixedpars=fixedpars,
                            dfns=dfns, aux=aux, cl=cl, sr.control = sr.control, ... = ...)

  })

  ## check on arg-types/values
  if (is.null(names(list_of_list_of_args)) || any( names(list_of_list_of_args) == "" ) )
    stop("Churn-types must be a named list.",call. = FALSE)
  stopifnot(is.data.frame(data))
  method <- match.arg(arg = method, choices = c('flexsurvreg','flexsurvspline'))

  ## type of censoring
  if (!is.null(time_lb_col_name) & !is.null(time_start_col_name) & is.element(method, c('flexsurvreg', 'flexsurvspline')) )
    stop(call. = FALSE,
         "\nYou've set a column both for the time lower-bound and for time-start, indicating data with both interval-censoring and truncation. ",
         "\nHowever, 'flexsurvreg' and 'flexsurvspline' currently only support one or the other of these, not both.")
  censoring_types <- c()
  if (is.null(time_lb_col_name)) time_lb_col_name <- time_col_name
  else censoring_types <- c(censoring_types, 'interval')
  if (!is.null(time_start_col_name)) censoring_types <- c(censoring_types, 'truncation')
  if (length(censoring_types)==0) censoring_types <- 'right'
  if ('interval' %in% censoring_types & method == 'flexsurvspline')
    stop(call. = FALSE, "Method 'flexsurvspline' does not currently support interval-censoring.")

  ## check on event-column:
  if (is.factor(data[[event_col_name]])) {
    censor_level <- levels(data[[event_col_name]])[[1]]
    if (censor_level %in% names(list_of_list_of_args))
      stop(call. = FALSE,
           "The first factor-level of the event column should correspond to the level indicating censoring; ",
           "however, ", censor_level, " was found in the names of `list_of_list_of_args`, which should only consist ",
           "of event-types. Please re-order this factor (e.g., using the `levels` function).")
  } else {
    censor_level <- setdiff(data[[event_col_name]], names(list_of_list_of_args))
    if (length(censor_level) != 1) {
      message("The levels of '",event_col_name,"' are: ", paste(unique(data[[event_col_name]]), collapse=", "), ".")
      message("However, the names of `list_of_list_of_args` are: ", paste(names(list_of_list_of_args), collapse=", "), ".")
      message("It should be the case that exactly one level in the event-column is not in the `list_of_list_of_args`-- this is the level indicating censoring.")
      if (length(censor_level)==0)
        message("If you want to allow for a censoring-level, but indicate that this particular dataset doesn't happen to have any censoring, ",
                "please set '", event_col_name, "' to a factor, and make the first level of this factor be the one indicating censoring.")
      stop(call. = FALSE, "Problem with the event-colum. See above.")
    }
  }
  missing_levels <- setdiff(unique(data[[event_col_name]]), c(censor_level, names(list_of_list_of_args)))
  if (length(missing_levels)>0)
    warning(call. = FALSE, immediate. = TRUE,
            "The following levels of ",event_col_name ," were not used in `list_of_list_of_args`: ",
            paste0(missing_levels, collapse = ", "))

  list_of_fits <- purrr::map2(
    .x = list_of_list_of_args,
    .y = names(list_of_list_of_args),
    function(args, this_event_name) {

      event_char <- paste0("(`", event_col_name, "`", "=='", this_event_name, "')")
      events_vec <- lazyeval::lazy_eval(as.lazy(event_char), data = data)
      if (all(!events_vec, na.rm = TRUE))
        warning(call. = FALSE,
                immediate. = TRUE,
                "There wasn't a single event of type '",this_event_name,"'.")

      if (method %in% c('flexsurvreg','flexsurvspline')) {
        if (censoring_types == 'interval')
          update_char <- paste0("tidysurv::interval_censor_helper(time= ", time_col_name, ", time_lb=", time_lb_col_name, ", event= ", event_char, ") ~ .")
        else if (censoring_types == 'truncation')
          update_char <- paste0("Surv(time=",time_start_col_name,", time2=", time_col_name,", event = ", event_char, ") ~ .")
        else
          update_char <- paste0("Surv(time=",time_col_name,", event = ", event_char, ") ~ .")
      }

      stopifnot( 'formula' %in% names(args) )
      args$formula <- update(args$formula, stats::as.formula(update_char))
      args$data <- data
      fit <- do.call( ifelse(memoize, yes = paste0("mem_",method), no = method) , args = args)
      fit$call$data <- NULL
      fit
    })

  class(list_of_fits) <- c('cr_survreg',class(list_of_fits))
  attr(list_of_fits,'cr_survreg') <- list(time_col_name = time_col_name,
                                          event_col_name = event_col_name,
                                          data = data,
                                          time_lb_col_name = time_lb_col_name,
                                          time_start_col_name = time_start_col_name,
                                          censor_level = censor_level,
                                          method = method,
                                          memoize = memoize)
  list_of_fits
}

#' @export
print.cr_survreg <- function(x, ...) {
  purrr::walk2(.x = x, .y = names(x), .f = function(ob,name) {
    cat("\n", name, "\n======\n", sep="")
    print(ob)
  })
}
#' @export
summary.cr_survreg <- function(object, ...) {
  print(object)
}



#' Create a tidy-surv object for an object of class \code{cr_survreg}
#'
#' @param object An object of class \code{cr_survreg}
#' @param newdata Optional. A dataframe.
#' @param group_vars A character-vector specifying column(s) to group-by when creating the
#'   survival-curves.
#' @param other_vars_fun All covariates not in \code{group_vars} will have this function applied to
#'   them. Should be a function that takes a numeric vector and returns a numeric vector with a single unique value.
#'   Defaults to `mean` with na-removal.
#' @param ... Ignored.
#'
#' @return An object of class \code{cr_survreg}
#' @export
#'
tidysurv.cr_survreg <- function(object, newdata = NULL, group_vars = NULL,
                                time_period = NULL,
                                other_vars_fun = function(x) mean(x, na.rm=TRUE),
                                ...) {
  if (is.null(newdata))
    newdata <- attr(object, "cr_survreg")$data

  terms_chars <- unique(purrr::flatten_chr(purrr::map(.x = object, .f = ~attr(delete.response(terms(.x)), "term.labels"))))
  df_mf <- model.frame(reformulate(terms_chars), data = newdata, na.action = na.pass)

  df_mapping <- purrr::map_df(object, ~get_terms_mappings(terms(.x), data = newdata), .id = 'sub_model')

  df_all_covs <- newdata[, unique(c(group_vars, df_mapping$variable_name)), drop = FALSE]

  # add any factors:
  df_to_add <- df_mapping %>% filter(term_class=='factor', !is.element(variable_name, group_vars))
  if (nrow(df_to_add)>0)
    for (i in seq_len(nrow(df_to_add))) {
      this_add <- df_to_add$column_name_in_mf[[i]]
      group_vars <- c(group_vars, this_add)
      df_all_covs[[this_add]] <- df_mf[[this_add]]
      message("Added '", df_to_add$column_name_in_mf[[i]], "' to the `group_vars` (factor).")
    }

  collapse_vars <- setdiff(colnames(df_all_covs), group_vars)
  newdata$..rownum <- 1:nrow(newdata)
  df_all_covs$..rownum <- 1:nrow(df_all_covs)
  if (length(group_vars) > 0)
    df_all_covs <- dplyr::group_by_(.data = df_all_covs, .dots = paste0("`", group_vars, "`"))

  if (length(collapse_vars) > 0)
    df_all_covs <- dplyr::mutate_at(.tbl = df_all_covs, .cols = vars(one_of(collapse_vars)), .funs = other_vars_fun)

  df_all_covs <- dplyr::ungroup(df_all_covs)

  ## get surv columns ---
  event_col_name <- attr(object,'cr_survreg')$event_col_name
  if (!is.element(event_col_name, colnames(newdata)))
    stop(call. = FALSE, "The column '", event_col_name, "' is not in `newdata`.")
  time_col_name <- attr(object,'cr_survreg')$time_col_name
  if (!all(df_all_covs$..rownum == newdata$..rownum))
    stop(call. = FALSE, "Please report this error to the package maintainer.") # i'm assuming dplyr will never rearrange due to grouping...
  df_all_covs$..rownum <- NULL
  df_all <- dplyr::bind_cols( newdata[,c(event_col_name,time_col_name),drop=FALSE], df_all_covs )

  # Create tidysurv Object ---
  censor_level <- attr(object,'cr_survreg')$censor_level
  df_all[[event_col_name]] <- factor(df_all[[event_col_name]], levels = c(censor_level, names(object)))
  form_lhs <- as.formula(paste0("Surv(time = `", time_col_name, "`, event = `", event_col_name, "`) ~ ."))

  if (length(group_vars)>0)
    df_tidysurv <- tidysurv(update(form_lhs, reformulate(paste0("`",group_vars,"`"))), data = df_all, time_period = time_period)
  else
    df_tidysurv <- tidysurv(update(form_lhs, .~1), data = df_all, time_period = time_period)

  # Get Predictions ---
  times <- sort(unique(df_tidysurv$time))
  df_small <- dplyr::distinct(df_all_covs)
  df_small$..rowname <- as.character(1:nrow(df_small))
  if (any(is.na(df_small)))
    stop(call. = FALSE, "NAs currently not allowed for this function. Please report this issue to the package maintainer.")
  list_of_fitted_dfs <- purrr::map2(
    .x = object,
    .y = names(object),
    .f = function(fit, event_name) {
      df_expanded <- tidyr::crossing(df_small, time = times)
      df_expanded[[event_name]] <-
        predict(fit, newdata = df_expanded,  type='survival', times = df_expanded$time)
      dplyr::as_data_frame(df_expanded[,c('time',event_name,'..rowname')])
    })
  df_covs_with_surv <- dplyr::left_join(df_small,
                                         purrr::reduce(list_of_fitted_dfs,dplyr::left_join, by=c('time','..rowname')),
                                         by = '..rowname')

  df_covs_with_surv$surv <- purrr::reduce(.x = df_covs_with_surv[names(object)], .f = `*`)
  df_covs_with_inc <- df_covs_with_surv %>%
    split(.$..rowname) %>%
    purrr::map_df(.f = function(df_chunk) {
      df_chunk[names(object)] <- purrr::map(
        .x = names(object),
        .f = function(event_name) {
          rate <- 1 - df_chunk[[event_name]]/lag(df_chunk[[event_name]], default = 1)
          cumsum(rate*lag(df_chunk$surv, default = 1))
        })
      return(df_chunk)
    })

  df_covs_with_haz <- dplyr::select(.data = df_covs_with_surv, -surv, -one_of(names(object)))
  df_covs_with_haz <-
    bind_cols(df_covs_with_haz,
              as.data.frame(predict(object, newdata = df_covs_with_haz, times = df_covs_with_haz$time, type = 'hazard')))
  df_covs_with_haz_tidy <- tidyr::gather_(data = df_covs_with_haz[,c('time',names(object),group_vars),drop=FALSE],
                                          key_col = '.surv_state',
                                          value_col = 'fitted_hazard',
                                          gather_cols = names(object))


  df_covs_with_inc <- df_covs_with_inc[,c('time',names(object),group_vars),drop=FALSE]
  df_fitted_final <- dplyr::distinct(tidyr::gather_(data = df_covs_with_inc,
                                                    key_col = '.surv_state',
                                                    value_col = 'fitted',
                                                    gather_cols = names(object)))
  df_fitted_final <- left_join(df_fitted_final, df_covs_with_haz_tidy, by = unname(c('time',group_vars,'.surv_state')))
  df_fitted_final$.surv_state <- factor(df_fitted_final$.surv_state, levels = levels(df_tidysurv$.surv_state))
  df_out <- dplyr::left_join(df_tidysurv, df_fitted_final, by = unname(c('time',group_vars,'.surv_state')))

  attr(df_out, 'tidysurv') <- list(strata_names = group_vars, method = 'cr_survreg', time_period = time_period)
  df_out

}

#' @export
predict.cr_survreg <- function(object, newdata = NULL, times, type = 'hazard', start = NULL, ...) {
  if (is.null(newdata))
    newdata <- attr(object,'cr_survreg')$data
  outl <- purrr::map(object, ~do.call(predict,
                                      args = c(
                                        list(.x, newdata = newdata, type = type, times = times, start = start),
                                        list(...))
  ))
  as.matrix(as.data.frame(outl))
}

#' Plot coefficients
#'
#' @export
plot_coefs <- function(object, ...) {
  UseMethod('plot_coefs')
}

#' @describeIn plot_coefs
#' Plot coefficients of `flexsurvreg` model
#'
#' @param object An object of type \code{cr_survreg}
#' @param ... Ignored
#' @return A ggplot object
#'
#' @export
plot_coefs.cr_survreg <- function(object, ...) {
  purrr::map(object, plot_coefs)
}

#' Create Surv object of type 'interval' given upper and lower bounds of event-times.
#'
#' @param time The time-column
#' @param time_lb The lower-bound for interval censoring.
#' @param event Event indicator (numeric or logical only)
#'
#' @return A \code{Surv} object with interval-censoring.
#' @export
interval_censor_helper <- function(time, time_lb, event) {
  args <- list()
  args$time <- ifelse(event==1, time_lb, time)
  args$time <- ifelse(args$time==0, .Machine$double.eps*10, args$time) # flexsurvspline complains about exact zeros
  args$time2 <- time
  args$event <- ifelse(time > time_lb, event*3, event)
  args$type <- "interval"
  do.call(survival::Surv, args)
}

#' Tidy an object of class `survfit`
#'
#' This function is/should be found in the `broom` package, but occasionally I want an immediate bug-
#' fix, and so I place it here.
#'
#' @param x An object of class `survfit`
#' @param ... Ignored.
#'
#' @return A tidy dataframe
#' @export
tidy.survfit <- function (x, ...)
{
  if (inherits(x, "survfitms")) {
    ret <- data.frame(time = x$time, n.risk = x$n.risk[,nrow(x$transitions),drop=TRUE],
                      n.event = c(x$n.event), n.censor = c(x$n.censor),
                      estimate = c(x$pstate), std.error = c(x$std.err),
                      conf.high = c(x$upper), conf.low = c(x$lower), state = rep(x$states,
                                                                                 each = nrow(x$pstate)))
    ret <- ret[ret$state != "", ]
  }
  else {
    ret <- data.frame(time = x$time, n.risk = x$n.risk, n.event = x$n.event,
                      n.censor = x$n.censor, estimate = x$surv, std.error = x$std.err,
                      conf.high = x$upper, conf.low = x$lower)
  }
  if (!is.null(x$strata)) {
    ret$strata <- rep(names(x$strata), x$strata)
  }
  ret
}

#' Get mapping from variables/column to formula terms
#'
#' @param terms A 'terms' object
#' @param data A dataframe
#'
#' @return A tidy dataframe with mapping.
#' @export
get_terms_mappings <- function(terms, data) {
  term_ob <- delete.response(terms)

  mapping_mat <- attr(term_ob, 'factors')
  mf_col_names_with_backticks <- rownames(mapping_mat)
  original_col_names <- purrr::map_chr(.x = mf_col_names_with_backticks, ~all.vars(lazyeval::as.lazy(.x)$expr))

  out <- data_frame(variable_name = original_col_names,
                    terms_with_variable = mf_col_names_with_backticks,
                    column_name_in_mf = stringr::str_replace(string = mf_col_names_with_backticks, pattern = "`", replacement = ""))

  model_frame <- model.frame(term_ob,data=data)
  out$term_class <- purrr::map_chr(.x = out$column_name_in_mf, .f = ~paste0(class(model_frame[[.x]]), collapse = "; "))
  out
}
