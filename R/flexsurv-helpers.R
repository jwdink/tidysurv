#' @export
terms.flexsurvreg <- function(x, ...) {
  stats::terms(x$concat.formula)
}

#' @describeIn plot_coefs
#' Plot coefficients of `flexsurvreg` model
#' @export
plot_coefs.flexsurvreg <- function(object, ...) {
  df_terms <-
    object$res.t[object$covpars,,drop=FALSE] %>%
    as.data.frame() %>% tibble::rownames_to_column('term') %>% as_data_frame()
  df_terms$parameter <- stringr::str_extract(string = df_terms$term, pattern = paste0(collapse="|", paste0(object$dlist$pars, "\\(")) )
  df_terms$parameter <- stringr::str_replace(string = df_terms$parameter, pattern = "\\(", replacement = "")
  df_terms$parameter <- coalesce(df_terms$parameter, object$dlist$location)

  df_terms %>%
    ggplot(aes(x = term, y = est)) +
    geom_pointrange(mapping = aes(ymin = `L95%`, ymax=`U95%`)) +
    facet_wrap(~parameter, scales = 'free') +
    geom_hline(yintercept = 0) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#' Form model-matrix
#'
#' Helper for \code{predict.flexsurvreg}.
#'
#' @param object Object of type \code{flexsurvreg}.
#' @param newdata A dataframe
#'
#' @return A model.matrix.
form.model.matrix <- function (object, newdata)
{
  mfo <- model.frame(object)
  covnames <- attr(mfo, "covnames")
  missing.covs <- unique(covnames[!covnames %in% names(newdata)])
  if (length(missing.covs) > 0) {
    missing.covs <- sprintf("\"%s\"", missing.covs)
    plural <- if (length(missing.covs) > 1)
      "s"
    else ""
    stop(sprintf("Value%s of covariate%s ", plural, plural),
         paste(missing.covs, collapse = ", "), " not supplied in \"newdata\"")
  }
  tt <- attr(mfo, "terms")
  Terms <- delete.response(tt)
  mf <- model.frame(Terms, newdata,na.action = na.omit,
                    xlev = .getXlevels(tt,mfo))
  if (!is.null(cl <- attr(Terms, "dataClasses")))
    .checkMFClasses(cl, mf)
  forms <- object$all.formulae
  mml <- vector(mode = "list", length = length(object$dlist$pars))
  names(mml) <- names(forms)
  forms[[1]] <- delete.response(terms(forms[[1]]))
  for (i in names(forms)) {
    mml[[i]] <- model.matrix(forms[[i]], mf)
  }

  ## compress model matrices
  cbind.drop.intercept <- function(...) do.call("cbind", lapply(list(...),
                                                                function(x) x[, -1, drop = FALSE]))
  X <- do.call("cbind.drop.intercept", mml)
  loc.cnames <- colnames(mml[[1]])[-1]
  anc.cnames <- unlist(mapply(function(x, y) sprintf("%s(%s)",
                                                     x, y), names(mml[-1]), lapply(mml[-1], function(x) colnames(x)[-1])))
  cnames <- c(loc.cnames, anc.cnames)
  colnames(X) <- cnames

  attr(X, "newdata") <- mf
  X
}

#' Add covs
#'
#' Helper for \code{predict.flexsurvreg}.
#'
#' @param x Object of type \code{flexsurvreg}.
#' @param pars Matrix of parameters
#' @param beta Matrix of betas
#' @param X Model-frame
#' @param transform Apply the inverse-transform associated with this family?
#'
#' @return Matrix of basepars
add.covs <- function (x, pars, beta, X, transform = FALSE)
{
  nres <- nrow(X)
  if (!is.matrix(pars))
    pars <- matrix(pars, nrow = nres, ncol = length(pars),
                   byrow = TRUE)
  if (!is.matrix(beta))
    beta <- matrix(beta, nrow = 1)
  for (j in seq(along = x$dlist$pars)) {
    covinds <- x$mx[[x$dlist$pars[j]]]
    if (length(covinds) > 0) {
      pars[, j] <- pars[, j] + beta[, covinds] %*% t(X[,covinds, drop = FALSE])
    }
    if (!transform)
      pars[, j] <- x$dlist$inv.transforms[[j]](pars[, j])
  }
  colnames(pars) <- x$dlist$pars
  pars
}

#' Predict method for `flexsurvreg` models
#'
#' This function returns predictions from a \code{flexsurvreg} object. This can be a
#' convenient/faster alternative to \code{summary.flexsurvreg}, when the number of observations with
#' unique combinations of covariates is large.
#'
#' @param object Object of type \code{flexsurvreg}.
#' @param newdata A dataframe to apply predictions to.
#' @param times A numeric vector of times for which to return predictions, with length equal to
#'   \code{nrow(newdata)}
#' @param type What type of prediction? Options are "survival", "cumhaz", and "hazard".
#' @param start Optional. A numeric vector of start (truncation) times.
#' @param ... Ignored.
#'
#' @return A numeric vector of predictions
#' @export
predict.flexsurvreg <- function(object, newdata, times, type = 'survival', start = NULL, ...) {

  if (length(times)==1)
    times <- rep(x = times, times = nrow(newdata))
  stopifnot( length(times) == nrow(newdata) )

  x <- object
  dat <- x$data
  Xraw <- model.frame(x)[, unique(attr(model.frame(x), "covnames.orig")), drop = FALSE]
  type <- match.arg(type, c("survival", "cumhaz", "hazard"))
  X <- form.model.matrix(object, as.data.frame(newdata))
  na_omit_idx <- attr(attr(X,'newdata'),'na.action')
  if (is.null(na_omit_idx)) na_omit_idx <- -(1:nrow(newdata))

  t <- times[-na_omit_idx]
  omit_idx <- which(t <= 0)
  t[omit_idx] <- median(t, na.rm = TRUE) # just a filler value

  if (is.null(start)) {
    start <- numeric(length(t))
  } else {
    if (length(start)==1)
      start <- rep(x = start, times = nrow(newdata))
    stopifnot( length(start) == nrow(newdata) )
    start <- start[-na_omit_idx]
  }


  fn <- switch(type, survival = function(t, start, ...) {
    ret <- (1 - x$dfns$p(t, ...))/(1 - x$dfns$p(start, ...))
    ret[t < start] <- 1
    ret
  }, hazard = function(t, start, ...) {
    ret <- x$dfns$h(t, ...) * (1 - x$dfns$p(start, ...))
    ret[t < start] <- 0
    ret
  }, cumhaz = function(t, start, ...) {
    ret <- x$dfns$H(t, ...) - x$dfns$H(start, ...)
    ret[t < start] <- 0
    ret
  })

  ##
  summfn2 <- fn
  args <- c(alist(t = , start = ), formals(fn))
  formals(summfn2) <- args[!duplicated(names(args))]
  body(summfn2) <- body(fn)
  fn <- summfn2

  ##
  fncall <- list(t, start)
  beta <- if (x$ncovs == 0) 0 else x$res[x$covpars, "est"]
  if ( (x$ncovs > 0) && (ncol(X) != length(beta)) ) {
    isare <- if (length(beta) == 1)
      "is"
    else "are"
    plural <- if (ncol(X) == 1)
      ""
    else "s"
    pluralc <- if (length(beta) == 1)
      ""
    else "s"
    stop("Supplied X has ", ncol(X), " column", plural,
         " but there ", isare, " ", length(beta), " covariate effect",
         pluralc)
  }

  dlist <- x$dlist

  basepars_mat <- add.covs(x, x$res.t[dlist$pars, "est"],
                           beta, X, transform=FALSE)
  basepars <- as.list(as.data.frame(basepars_mat))
  fncall[dlist$pars] <- basepars[dlist$pars]
  fncall[names(x$aux)] <- x$aux
  y <- do.call(fn, fncall)
  if (type == 'survival') y[omit_idx] <- 1
  else y[omit_idx] <- 0

  # fill for NAs
  out <- numeric(nrow(newdata))
  out[-na_omit_idx] <- y
  if (!all(na_omit_idx>0))
    out[na_omit_idx] <- NA_real_

  return(out)


}

#
# mem_flexsurvspline <- memoise::memoise(function(formula, data, weights=NULL,bhazard=NULL,subset=NULL,k=0,knots=NULL,bknots=NULL,scale='hazard',timescale='log',...) {
#   if ( (!is.null(weights)) | (!is.null(bhazard)) | (!is.null(subset)) )
#     stop("Weights, bhazard, and subset are non supported if `memoize= TRUE`.")
#   flexsurv::flexsurvspline(formula = formula, data = data, k = k, knots = knots, bknots = bknots, scale = scale, timescale = timescale, ... = ...)
# })
#

