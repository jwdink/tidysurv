% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidysurv-package.R
\name{cr_survreg}
\alias{cr_survreg}
\title{Easy interface for parametric competing-risks survival-regression}
\usage{
cr_survreg(time_col_name, event_col_name, data, list_of_list_of_args,
  time_lb_col_name = NULL, time_start_col_name = NULL,
  method = "flexsurvreg", memoize = TRUE)
}
\arguments{
\item{time_col_name}{A character indicating the column-name of the 'time' column.}

\item{event_col_name}{A character indicating the column-name for events. The values in this
column (either factor or character), should match the names of the \code{list_of_list_of_args},
except for the value indicating censoring.}

\item{data}{A dataframe.}

\item{list_of_list_of_args}{A list, whose names correspond to each possible type of event. These
should correspond to the values in the \code{event_col_name} column. The elements of this list
are themselves lists— arguments to be passed to the survival-regression modelling function (as
specified by \code{method}— currently only \code{flexsurvreg} is supported).}

\item{time_lb_col_name}{Optional. A character indicating the column-name of the time
'lower-bound'. See the vignette for details.}

\item{time_start_col_name}{Optional. A character indicating the column-name of the start
(truncation) times. See the vignette for details.}

\item{method}{A character string naming the function to be called for survival-regression
modelling. Currently only supports \code{flexsurvreg}.}

\item{memoize}{When this function runs, it calls the fitting method (e.g., \code{flexsurvreg}),
once per event-type. If you're only tweaking the model for one event-type, you might not want to
waste time re-fitting the model for all of the other event types. When \code{memoize} is TRUE (the
default), calls to the fitting method are cached: if called with the exact same arguments, the
code isn't re-run. This means that you can tweak the model for one event-type at a time, without
having to worry about re-fitting the other models on each function call.}
}
\value{
An object of type \code{cr_survreg}, with plot and summary methods.
}
\description{
Easy interface for parametric competing-risks survival-regression
}
