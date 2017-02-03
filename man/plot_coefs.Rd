% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flexsurv-helpers.R, R/tidysurv-package.R
\name{plot_coefs.flexsurvreg}
\alias{plot_coefs}
\alias{plot_coefs.cr_survreg}
\alias{plot_coefs.flexsurvreg}
\title{#' Title
#'
#' @param object
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
msfit <- function(object, ...) {
  UseMethod('msfit')
}}
\usage{
\method{plot_coefs}{flexsurvreg}(object, ...)

plot_coefs(object, ...)

\method{plot_coefs}{cr_survreg}(object, ...)
}
\arguments{
\item{object}{An object of type \code{flexsurvreg}}

\item{...}{Ignored}

\item{object}{An object of type \code{cr_survreg}}

\item{...}{Ignored}
}
\value{
A ggplot object

A ggplot object
}
\description{
#' Title
#'
#' @param object
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
msfit <- function(object, ...) {
  UseMethod('msfit')
}

Plot coefficients
}
\section{Methods (by class)}{
\itemize{
\item \code{flexsurvreg}: Plot coefficients of `flexsurvreg` model

\item \code{cr_survreg}: Plot coefficients of `flexsurvreg` model
}}
