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
#' @importFrom matrixStats rowCollapse
#' @importFrom glue glue
#' @importFrom tibble rownames_to_column
#' @importFrom broom tidy
#' @importFrom purrr map map2 pmap map_dbl map_chr map_lgl
#'             map2_chr map2_dbl map2_lgl transpose flatten_chr flatten_dbl flatten_lgl
#'             walk walk2 map_df map2_df
#' @importFrom graphics plot
#' @importFrom stats sd terms update integrate delete.response
#'             as.formula coef fitted glm lm median na.omit poly
#'             predict var quantile model.response model.frame
#'             na.pass na.exclude na.fail model.matrix model.weights
#'             .getXlevels .checkMFClasses reformulate binom.test formula
#'             setNames
#' @importFrom grDevices dev.off pdf
#'
#' @docType package
#' @name tidysurv
NULL
#> NULL

#' Take a dataframe, and convert its time column(s) to binned versions.
#'
#' @param data A data.frame
#' @param time_period Number for bin-unit
#' @param formula_lhs The expression on the left-hand-side of the formula, which generates the `Surv` object
#'
#' @return The data, with the relevant time-columns now binned.
#' @export
convert_time_cols_to_binned <- function(data, time_period, formula_lhs) {
  stopifnot(is.data.frame(data))
  stopifnot(length(time_period)==1)

  response_object <- eval(formula_lhs, envir = data)
  if (class(response_object)[1]=="Survint")
    response_object <- with(as.data.frame(response_object), Surv(time = start, time2 = end, event = event))
  bin_times <- function(time_vec, time_period, start = time_period) {
    time_grid <- seq(from = start, to = max(time_vec,na.rm=TRUE), by = time_period)
    time_mat <- matrix(data = rep(time_grid, each = length(time_vec)), nrow = length(time_vec))
    matrixStats::rowCollapse(x = time_mat, idxs = apply(X = abs(time_vec - time_mat), FUN = which.min, MARGIN = 1))
  }
  if ( attr(response_object, 'type') %in% c('counting','mcounting') ) {
    df_response <- as.data.frame(as.matrix(response_object))
    min_diff <- min(df_response$stop - df_response$start)
    if (time_period > min_diff)
      stop(call. = FALSE,
           "It looks like you're using 'counting' format for your DV (i.e., you have time-dependent covariates).\n",
           "In this case, `time_period` cannot be less than the smallest time-interval (min(stop-start)), which is ", min_diff)
    time_col_name <- as.character(formula_lhs$time2)
    time_start_col_name <- as.character(formula_lhs$time)
    if ( !is.element(time_col_name, colnames(data)) | !is.element(time_start_col_name, colnames(data)) )
      stop(call. = FALSE,
           "If using `time_period`, your `Surv` response object should only reference column-names (not expressions including column-names).")
    data[[time_col_name]] <- bin_times(data[[time_col_name]], time_period)
    data[[time_start_col_name]] <- bin_times(data[[time_start_col_name]], time_period, start = min(data[[time_start_col_name]]))
  } else {
    time_col_name <- as.character(formula_lhs$time)
    if (!is.element(time_col_name, colnames(data)))
      stop(call. = FALSE,
           "If using `time_period`, your `Surv` response object should only reference column-names (not expressions including column-names).")
    data[[time_col_name]] <- bin_times(data[[time_col_name]], time_period)
  }
  return(data)
}

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


  formula_lhs <- lazyeval::call_standardise(object[[2]])

  if (!is.null(time_period))
    data <- convert_time_cols_to_binned(data = data, time_period = time_period, formula_lhs = formula_lhs)
  response_object <- eval(formula_lhs, envir = data)
  if (class(response_object)[1]=="Survint")
    response_object <- with(as.data.frame(response_object), Surv(time = start, time2 = end, event = event))

  # ## create model-frame, mapping:
  model_frame <- model.frame(formula = object, data = data, na.action=na.pass)
  model_frame <- model_frame[,-1,drop=FALSE]

  # call survfit on pre-stratified data:
  if (ncol(model_frame)>0) {
    df_mf_distinct <- dplyr::distinct(model_frame)
    strata_names <- colnames(df_mf_distinct)
    df_mf_distinct <- rownames_to_column(df_mf_distinct, var = '..strata')
    if (nrow(df_mf_distinct) > max_num_levels)
      stop(call. = FALSE, "The right-hand-side of the formula defines more than ",
           max_num_levels,
           " unique levels. Increase `max_num_levels` if this was intentional.")
    df_mf_w_strata <- left_join(model_frame, df_mf_distinct, by = colnames(model_frame))
    sfit <- survival::survfit(formula = update(response_object~1, .~..strata), data = df_mf_w_strata,
                    na.action = na.action, ... = ...)
  } else {
    sfit <- survival::survfit(formula = response_object~1, na.action = na.action, ... = ...)
    strata_names <- NULL
  }

  ## convert strata -> covariates:
  df_tidy <- tidy.survfit(sfit)
  if (inherits(sfit, "survfitms"))
    df_tidy <- dplyr::rename(.data = df_tidy, .surv_state = state)
  if ( ncol(model_frame)>0 )  {
    df_tidy$..strata <- stringr::str_match(string = df_tidy$strata, pattern = "\\.\\.strata=(.*)")[,2]
    out <- dplyr::left_join(x = df_tidy, y = df_mf_distinct, by = '..strata')
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
      g <- g + geom_line(mapping = aes(y = fitted_rate), size=1.2,alpha=.75)
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
#'   are themselves lists- arguments to be passed to the survival-regression modelling function (as
#'   specified by \code{method}- currently only \code{flexsurvreg} is supported).
#' @param time_lb_col_name Optional. A character indicating the column-name of the time
#'   'lower-bound'. See the vignette for details.
#' @param time_dependent_config If your data is in a 'counting process' format (i.e., there are
#'   multiple rows per 'person', with a 'start' column specifying right-truncation), you should
#'   supply this. A list with three entries: `time_start_col_name`, `id_col_name`, and
#'   `time_dependent_col_names`.
#' @param method A character string naming the function to be called for survival-regression
#'   modelling. Currently only supports \code{flexsurvreg}.
#'
#' @return An object of type \code{cr_survreg}, with plot and summary methods.
#' @export
#'
cr_survreg <- function(time_col_name,
                       event_col_name,
                       data,
                       list_of_list_of_args,
                       time_lb_col_name = NULL,
                       time_dependent_config = list(time_start_col_name = NULL, id_col_name = NULL, time_dependent_col_names = NULL),
                       method = 'flexsurvreg') {

  ## check on arg-types/values
  if (is.null(names(list_of_list_of_args)) || any( names(list_of_list_of_args) == "" ) )
    stop("Churn-types must be a named list.",call. = FALSE)
  stopifnot(is.data.frame(data))
  method <- match.arg(arg = method, choices = c('flexsurvreg','flexsurvspline','survreg_map'))

  stopifnot( c('time_start_col_name', 'id_col_name', 'time_dependent_col_names') %in% names(time_dependent_config))
  time_start_col_name <- time_dependent_config$time_start_col_name

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
      event_char <- paste0("as.integer(`", event_col_name, "`", "=='", this_event_name, "')")
      events_vec <- lazyeval::lazy_eval(as.lazy(event_char), data = data)
      if (all(!events_vec, na.rm = TRUE))
        warning(call. = FALSE,
                immediate. = TRUE,
                "There wasn't a single event of type '",this_event_name,"'.")

      if (method %in% c('flexsurvreg','flexsurvspline')) {
        if (censoring_types == 'interval')
          update_char <- glue::glue("interval_censor_helper(time={time_col_name}, time_lb = {time_lb_col_name}, event = {event_char}) ~ .")
        else if (censoring_types == 'truncation')
          update_char <- glue::glue("Surv(time={time_start_col_name}, time2={time_col_name}, event={event_char}) ~ .")
        else
          update_char <- glue::glue("Surv(time={time_col_name}, event={event_char}) ~ .")
      } else if (method == 'survreg_map') {
        the_call <- call('Survint',
                         end = parse(text = time_col_name)[[1]],
                         event = parse(text = event_char)[[1]])
        if (!is.null(time_start_col_name)) the_call$start <- parse(text = time_start_col_name)[[1]]
        if (!is.null(time_lb_col_name)) the_call$end_lb <- parse(text = time_lb_col_name)[[1]]
        update_char <- paste(paste0(deparse(the_call), collapse=""), "~.")
      }
      stopifnot( 'formula' %in% names(args) )
      args$formula <- update(args$formula, stats::as.formula(update_char))
      args$data <- data
      fit <- do.call( method , args = args)
      fit$call$data <- NULL
      fit
    })

  if (!is.null(time_dependent_config$time_dependent_col_names))
    walk(time_dependent_config$time_dependent_col_names,
         ~if (!is.numeric(data[[.x]])) stop(call. = FALSE, "Currently, only numeric time-dependent variables are supported."))

  class(list_of_fits) <- c('cr_survreg',class(list_of_fits))
  attr(list_of_fits,'cr_survreg') <- list(time_col_name = time_col_name,
                                          event_col_name = event_col_name,
                                          data = data,
                                          id_col_name = time_dependent_config$id_col_name,
                                          time_dependent_col_names = time_dependent_config$time_dependent_col_names,
                                          time_lb_col_name = time_lb_col_name,
                                          time_start_col_name = time_start_col_name,
                                          censor_level = censor_level)
  list_of_fits
}

#' Create a `cr_survreg` model from a custom list of models
#'
#' Instead of using a predetermined type of survival-regression model for the submodels in
#' cr_survreg, you can use any survival-regression model with a valid `predict` method (one that
#' takes \code{type='survival'} and \code{times}; see \code{?predict.flexsurvreg}).
#'
#' @param list_of_models  A list, whose names correspond to each possible type of event. These
#'   should correspond to the values in the \code{event_col_name} column. Each element is a
#'   survival-regression model that can predict survival probabilities.
#' @param data A data.frame
#' @param time_col_name A character indicating the column-name of the 'time' column.
#' @param event_col_name A character indicating the column-name for events. The values in this
#'   column (either factor or character), should match the names of the names of
#'   \code{list_of_models}, except for the value indicating censoring.
#' @param time_lb_col_name Optional. A character indicating the column-name of the time
#'   'lower-bound'.
#' @param time_dependent_config If your data is in a 'counting process' format (i.e., there are
#'   multiple rows per 'person', with a 'start' column specifying right-truncation), you should
#'   supply this. A list with three entries: `time_start_col_name`, `id_col_name`, and
#'   `time_dependent_col_names`.
#'
#' @return An object of type \code{cr_survreg}, with plot and summary methods.
#' @export
convert_to_cr_survreg <- function(list_of_models,
                                  data,
                                  time_col_name,
                                  event_col_name,
                                  time_lb_col_name = NULL,
                                  time_dependent_config = list(time_start_col_name = NULL, id_col_name = NULL, time_dependent_col_names = NULL)) {

  stopifnot( c('time_start_col_name', 'id_col_name', 'time_dependent_col_names') %in% names(time_dependent_config))
  time_start_col_name <- time_dependent_config$time_start_col_name

  ## check on event-column:
  if (is.factor(data[[event_col_name]])) {
    censor_level <- levels(data[[event_col_name]])[[1]]
    if (censor_level %in% names(list_of_models))
      stop(call. = FALSE,
           "The first factor-level of the event column should correspond to the level indicating censoring; ",
           "however, ", censor_level, " was found in the names of `list_of_list_of_args`, which should only consist ",
           "of event-types. Please re-order this factor (e.g., using the `levels` function).")
  } else {
    censor_level <- setdiff(data[[event_col_name]], names(list_of_models))
    if (length(censor_level) != 1) {
      message("The levels of '",event_col_name,"' are: ", paste(unique(data[[event_col_name]]), collapse=", "), ".")
      message("However, the names of `list_of_models` are: ", paste(names(list_of_models), collapse=", "), ".")
      message("It should be the case that exactly one level in the event-column is not in the `list_of_models`-- this is the level indicating censoring.")
      if (length(censor_level)==0)
        message("If you want to allow for a censoring-level, but indicate that this particular dataset doesn't happen to have any censoring, ",
                "please set '", event_col_name, "' to a factor, and make the first level of this factor be the one indicating censoring.")
      stop(call. = FALSE, "Problem with the event-colum. See above.")
    }
  }
  missing_levels <- setdiff(unique(data[[event_col_name]]), c(censor_level, names(list_of_models)))
  if (length(missing_levels)>0)
    warning(call. = FALSE, immediate. = TRUE,
            "The following levels of ",event_col_name ," were not used in `list_of_models`: ",
            paste0(missing_levels, collapse = ", "))


  if (!is.null(time_dependent_config$time_dependent_col_names))
    walk(time_dependent_config$time_dependent_col_names,
         ~if (!is.numeric(data[[.x]])) stop(call. = FALSE, "Currently, only numeric time-dependent variables are supported."))


  class(list_of_models) <- c('cr_survreg',class(list_of_models))
  attr(list_of_models,'cr_survreg') <- list(time_col_name = time_col_name,
                                          event_col_name = event_col_name,
                                          data = data,
                                          id_col_name = time_dependent_config$id_col_name,
                                          time_dependent_col_names = time_dependent_config$time_dependent_col_names,
                                          time_lb_col_name = time_lb_col_name,
                                          time_start_col_name = time_start_col_name,
                                          censor_level = censor_level)
  list_of_models

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


#' Merge formulae
#'
#' @param forms List of formulae
#' @param data Data.frame
#'
#' @return Merged formula
#' @export
merge_formulae <- function(forms, data) {
  term_objs <- purrr::map(forms, ~terms(.x, data=data))
  covnames <- unique(purrr::flatten_chr( purrr::map(term_objs, ~attr(.x,'term.labels')) ))
  cov_formula_char <- if (length(covnames)==0) "1" else paste0(collapse = " + ", covnames)
  out <- as.formula(paste0("~", cov_formula_char))
  environment(out) <- environment(forms[[1]])
  out
}

#' Create a tidy-surv object for an object of class \code{cr_survreg}
#'
#' @param object An object of class \code{cr_survreg}
#' @param newdata Optional. A dataframe.
#' @param group_vars A character-vector specifying column(s) to group-by when creating the
#'   survival-curves.
#' @param time_period This bins the time-values into groups of this width, which is needed for
#'   plotting rates and time-dependent variables.
#' @param id_col_name Character-string of id-column. Only needed if originally supplied value
#'   doesn't work for `newdata`
#' @param time_dependent_vars Character-strings for time-dependent covariates. Instead of taking the
#'   overall mean, takes the mean at each time-point.
#' @param ... Passed to tidysurv.formula, which in turn passes these args to \code{survfit}
#'
#' @return An object of class \code{cr_survreg}
#' @export
#'
tidysurv.cr_survreg <- function(object, newdata = NULL, group_vars = NULL,
                                time_period = NULL,
                                id_col_name = NULL, time_dependent_vars = NULL,
                                max_num_levels = 200,
                                ...) {

  if (is.null(newdata))
    newdata <- attr(object, "cr_survreg")$data
  newdata <- ungroup(newdata)

  ## create model-frame, mapping:
  forms <- map(map(map(object,terms), delete.response), formula)
  full_formula <- merge_formulae(forms, data = newdata)
  terms_mapper <- get_terms_mapper(formula = full_formula, data = newdata)
  model_frame <- model.frame(formula = full_formula, data = newdata, na.action=na.pass)
  numeric_lgl <- map_lgl(model_frame, is.numeric)
  from_mf_to_od <- terms_mapper(model_frame_cols = names(numeric_lgl))
  vars_to_add_to_group <- names(which(!numeric_lgl))

  ## identify non-numeric variables (can't be collapsed)
  # we have to do collapsing on the original data (not model-frame),
  # because the former is what `predict` takes. but we can identify
  # rows to group-by with the model-frame.
  df_covariates <- bind_cols(
    newdata[, unique(c(group_vars, flatten_chr(from_mf_to_od))), drop=FALSE],
    model_frame[,vars_to_add_to_group,drop=FALSE])

  group_vars <- unique(c(group_vars, vars_to_add_to_group))

  # deal with time-period:
  attrs <- attr(object, 'cr_survreg')
  if (!is.null(attrs$time_start_col_name))
    f_lhs_char <- glue::glue("Surv(time = {time_start_col_name}, time2 = {time_col_name}, event = {event_col_name})",
                             .envir = as.environment(attrs))
  else
    f_lhs_char <- glue::glue("Surv(time = {time_col_name}, event = {event_col_name})", .envir = as.environment(attrs))
  formula_lhs <- parse(text=f_lhs_char)[[1]]
  time_col_name <- attr(object,'cr_survreg')$time_col_name
  if (!is.factor(newdata[[attrs$event_col_name]]))
    newdata[[attrs$event_col_name]] <- relevel(as.factor(newdata[[attrs$event_col_name]]),
                                                           ref = attrs$censor_level)
  if (!is.null(time_period))
    newdata <- convert_time_cols_to_binned(data = newdata, time_period = time_period, formula_lhs = formula_lhs)

  if (is.null(time_dependent_vars))
    time_dependent_vars <- attrs$time_dependent_col_names
  else if (is.na(time_dependent_vars)) # this explicitly overrides them
    time_dependent_vars <- NULL
  else
    if (is.null(attrs$time_start_col_name))
      stop("You've specified `time_dependent_vars`, but this model isn't in 'counting process' format ",
           "(i.e., there was no 'time_start_col_name' found in the model-object).")

  ## group by group vars, collapse all others.
  collapse_vars <- setdiff(colnames(df_covariates), group_vars)
  collapse_vars <- setdiff(collapse_vars, time_dependent_vars)
  df_covariates <- rownames_to_column(df_covariates, var = '..row')
  if (is.null(id_col_name)) id_col_name <- attrs$id_col_name
  if (!is.null(id_col_name)) df_covariates[[id_col_name]] <- newdata[[id_col_name]]
  if (length(group_vars)>0)
    split_df_covs <- split(x = df_covariates[,c('..row',collapse_vars,id_col_name),drop=FALSE],
                           f = as.list(df_covariates[,group_vars,drop=FALSE]), drop=TRUE)
  else
    split_df_covs <- list(df_covariates[,c('..row',collapse_vars,id_col_name),drop=FALSE])


  if (length(collapse_vars)>0) {
    df_collapsed <- map_df(split_df_covs,
                           function(df_chunk) {
                             if (!is.null(id_col_name))
                               df_start <- df_chunk %>%
                                 group_by_(.dots = id_col_name) %>%
                                 summarize_at(.cols = collapse_vars, .funs = mean, na.rm=TRUE)
                             else
                               df_start <- df_chunk
                             df_means <- df_start %>%
                               summarize_at(.cols = collapse_vars, .funs = mean, na.rm=TRUE)
                             for (col in colnames(df_means))
                               df_chunk[[col]] <- df_means[[col]]
                             as_data_frame(df_chunk)
                           })
    df_collapsed <- left_join(x = df_covariates[,c('..row',group_vars),drop=FALSE],
                              y = df_collapsed,
                              by = '..row')
  } else {
    df_collapsed <- left_join(x = df_covariates[,c('..row',group_vars),drop=FALSE],
                              y = bind_rows(split_df_covs),
                              by = '..row')
  }
  all_times <- newdata[[time_col_name]]
  if (!is.null(attrs$time_start_col_name)) all_times <- c(all_times, newdata[[attrs$time_start_col_name]])
  if (is.null(time_period))
    unique_times <- sort(unique(all_times))
  else
    unique_times <- full_seq(x = all_times, period = time_period)

  ## if there are time-dependent covariates, we'll be using the averages
  # for each time-point, not overall averages. these will be calculated now,
  # then used in prediction in the next step
  if (!is.null(time_dependent_vars)) {
    # TO DO: degenerative behavior when a variable doesn't actually vary with time
    # might help you understand if this is really doing what it should.
    if (is.null(time_period))
      stop(call. = FALSE, "If there are time-dependent covariates, you must specify `time_period`.")
    df_td <- newdata[,c(time_dependent_vars,id_col_name,attrs$time_start_col_name),drop=FALSE]
    df_grid <- crossing(ID = unique(df_td[[id_col_name]]),TIME = unique_times)
    df_grid[[id_col_name]] <- df_grid$ID; df_grid[[attrs$time_start_col_name]] <- df_grid$TIME;
    df_grid$ID <- df_grid$TIME <- NULL
    df_td_full <- left_join(df_grid, df_td, by = c(id_col_name, attrs$time_start_col_name))
    df_max_times <- newdata %>%
      group_by_(.dots = id_col_name) %>%
      summarize_(.dots = list(..max = lazyeval::interp(~max(TIME), TIME = as.name(time_col_name))))
    df_td_full <- left_join(df_td_full, df_max_times, by = id_col_name)
    df_td_full <- df_td_full[df_td_full[[attrs$time_start_col_name]] <= df_td_full$..max, ]
    df_td_full <- arrange_(df_td_full, .dots = c(id_col_name, attrs$time_start_col_name))
    df_td_full <- df_td_full %>%
      group_by_(.dots = id_col_name) %>%
      fill_(fill_cols = time_dependent_vars, .direction = 'down') %>%
      ungroup()

    if (length(group_vars)>0) {
      if (any(group_vars %in% time_dependent_vars))
        stop(call. = FALSE, "None of the `group_vars` can be time-dependent variables.")
      df_id_group_vals <- distinct(df_covariates[,c(id_col_name, group_vars),drop=FALSE])
      df_td_full <- left_join(x = df_td_full, y = df_id_group_vals, by = id_col_name)
    }
    df_td_values <- df_td_full %>%
      group_by_(.dots = c(attrs$time_start_col_name, group_vars)) %>%
      summarize_at(.cols = time_dependent_vars, .funs = mean, na.rm=TRUE) %>%
      ungroup()
  }

  # Create tidysurv Object ---
  if (!is.null(attrs$time_start_col_name))
    df_collapsed_w_resp <- bind_cols(df_collapsed,
                                 newdata[,flatten_chr(attrs[c('time_col_name','event_col_name','time_start_col_name')]),drop=FALSE])
  else
    df_collapsed_w_resp <- bind_cols(df_collapsed,
                                 newdata[,flatten_chr(attrs[c('time_col_name','event_col_name')]),drop=FALSE])

  ts_form <- as.formula(paste(f_lhs_char, "~1"))
  if (length(group_vars)>0)
    ts_form <- update(ts_form, reformulate(paste0("`",group_vars,"`")))
    df_tidysurv <- tidysurv(ts_form, data = df_collapsed_w_resp, time_period = time_period, max_num_levels=max_num_levels,...=...)

  # Get Predictions ---
  if (!is.null(id_col_name)) df_collapsed[[id_col_name]] <- NULL
  if (ncol(df_collapsed)==1) {
    df_small <- data_frame(..row = '1')
  } else {
    df_small <- dplyr::distinct(select(df_collapsed, -`..row`))
    df_small <- rownames_to_column(df_small, var = "..row")
  }
  df_small <- na.exclude(df_small)

  df_expanded <- tidyr::crossing(df_small, time = unique_times)
  if (!is.null(time_dependent_vars))
    df_expanded <- left_join(df_expanded, df_td_values,
                             by = c(`time` = attrs$time_start_col_name, group_vars))
  if (length(group_vars)>0)
    df_expanded <- group_by_(.data = df_expanded, .dots = group_vars)
  df_expanded <- df_expanded %>%
    mutate(time_lagged = lag(time, default = 0)) %>%
    ungroup() %>%
    as_data_frame()
  list_of_fitted_dfs <- purrr::map2(
    .x = object,
    .y = names(object),
    .f = function(fit, event_name) {
      surv <- predict(fit, newdata = df_expanded,  type='survival', times = df_expanded$time)
      surv_lagged <- predict(fit, newdata = df_expanded,  type='survival', times = df_expanded$time_lagged)
      df_expanded[[event_name]] <- 1-surv/surv_lagged
      dplyr::as_data_frame(df_expanded[,c('time',event_name,'..row')])
    })

  df_covs_with_rate <- dplyr::left_join(df_small,
                                         purrr::reduce(list_of_fitted_dfs,dplyr::left_join, by=c('time','..row')),
                                         by = '..row')
  df_covs_with_rate$all_churn <- rowSums(df_covs_with_rate[,names(object),drop=FALSE], na.rm = FALSE)
  df_covs_with_rate <- df_covs_with_rate %>%
    group_by(..row) %>%
    mutate(surv = cumprod(1-all_churn), all_churn=NULL) %>%
    ungroup()
  df_covs_with_inc <- df_covs_with_rate %>%
    split(.$..row) %>%
    purrr::map_df(.f = function(df_chunk) {
      df_chunk[names(object)] <- purrr::map(
        .x = names(object),
        .f = function(event_name) {
          cumsum(df_chunk[[event_name]]*lag(df_chunk$surv, default = 1))
        })
      return(df_chunk)
    })

  df_covs_with_rate_tidy <- tidyr::gather_(data = df_covs_with_rate[,c('time',names(object),group_vars),drop=FALSE],
                                          key_col = '.surv_state',
                                          value_col = 'fitted_rate',
                                          gather_cols = names(object))


  df_covs_with_inc <- df_covs_with_inc[,c('time',names(object),group_vars),drop=FALSE]
  df_fitted_final <- dplyr::distinct(tidyr::gather_(data = df_covs_with_inc,
                                                    key_col = '.surv_state',
                                                    value_col = 'fitted',
                                                    gather_cols = names(object)))
  df_fitted_final <- left_join(df_fitted_final, df_covs_with_rate_tidy, by = unname(c('time',group_vars,'.surv_state')))
  df_fitted_final$.surv_state <- factor(df_fitted_final$.surv_state, levels = levels(df_tidysurv$.surv_state))
  df_out <- dplyr::left_join(df_tidysurv, df_fitted_final, by = unname(c('time',group_vars,'.surv_state')))

  attr(df_out, 'tidysurv') <- list(strata_names = group_vars, method = 'cr_survreg', time_period = time_period)
  df_out

}

#' @export
predict.cr_survreg <- function(object, newdata = NULL, times, type = 'survival', start = NULL, ...) {
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
#' @param object An object
#' @param ... Passed to methods
#'
#' @return A ggplot object
#'
#' @export
plot_coefs <- function(object, ...) {
  UseMethod('plot_coefs')
}

#' @describeIn plot_coefs
#' Plot coefficients of `cr_survreg` model
#' @export
plot_coefs.cr_survreg <- function(object, ...) {
  purrr::map(object, plot_coefs, ... = ...)
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

#' Get a function which maps from column names to model-terms and vice versa
#'
#' @param formula The formula to be passed to model.frame
#' @param data The data.frame
#' @param ... Arguments to be passed to \code{stats::model.matrix}
#'
#' @return A function that gives the model-term(s) for each original-column, or gives the
#'   original-column(s) for each model-term. You can also get original column-names from model-frame
#'   names (e.g., instead of specific levels of factor, just factor(variable) gets mapped to
#'   variable).
#' @export
get_terms_mapper <- function(formula, data, ...) {
  # remove response:
  if (length(formula)>2) formula[[2]] <- NULL
  model_frame <- model.frame(formula = formula, data, na.action=na.pass)
  if (ncol(model_frame)==0) {
    from_mm_to_od <- from_od_to_mm <- from_mf_to_od <- list()
  } else {
    model_matrix <- do.call(model.matrix, c(list(terms(model_frame), data=model_frame), list(...)))

    cols_in_mm <- colnames(model_matrix)
    fact_mat <- attr(terms(model_frame),'factors')
    cols_in_mf <- map(attr(model_matrix,'assign'),
                      function(assign_idx) row.names(fact_mat)[1==fact_mat[,assign_idx,drop=TRUE]])

    ##
    is_null_lgl <- map_lgl(cols_in_mf, is.null)
    cols_in_od <- vector(mode = 'list', length = length(is_null_lgl))
    if (any(!is_null_lgl)) {
      cols_in_od[!is_null_lgl] <- map(
        .x = cols_in_mf[!is_null_lgl],
        .f = function(vec_of_mm_cols) flatten_chr(map(vec_of_mm_cols, ~all.vars(parse(text = .x)[[1]]))))
      cols_in_od[!is_null_lgl] <- map(cols_in_od[!is_null_lgl], ~.x[.x%in%colnames(data)])
    }

    ##
    from_mf_to_od <- map(colnames(model_frame), ~all.vars(parse(text = .x)[[1]])) %>%
      map(~.x[.x%in%colnames(data)])
    names(from_mf_to_od) <- colnames(model_frame)

    ##
    from_mm_to_od <- setNames(nm = cols_in_mm, cols_in_od)
    from_od_to_mm <- data_frame(model_mat_col = cols_in_mm, original_col = cols_in_od) %>%
      unnest() %>%
      group_by(original_col) %>%
      do(model_mat_cols = .$model_mat_col) %>%
      ungroup() %>%
      deframe()
  }

  function(original_cols = NULL, model_mat_cols = NULL, model_frame_cols = NULL) {
    if (!is.null(model_mat_cols))
      from_mm_to_od[model_mat_cols]
    else if (!is.null(original_cols))
      from_od_to_mm[original_cols]
    else if (!is.null(model_frame_cols))
      from_mf_to_od[model_frame_cols]
  }
}

