% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_moeadps.R
\name{plot.moeadps}
\alias{plot.moeadps}
\title{plot.moeadps}
\usage{
\method{plot}{moeadps}(
  x,
  ...,
  useArchive = FALSE,
  feasible.only = TRUE,
  viol.threshold = 1e-06,
  nondominated.only = TRUE,
  plot.weights = FALSE,
  which.objectives = NULL,
  suppress.pause = FALSE,
  color.by.obj = 1
)
}
\arguments{
\item{x}{list object of class \emph{moead}
(generated by \code{\link[=moead]{moead()}})}

\item{...}{other parameters to be passed down to specific plotting functions (currently unused)}

\item{useArchive}{logical flag to use information from \code{x$Archive}.
Only used if x$Archive is not \code{NULL}.}

\item{feasible.only}{logical flag to use only feasible points in the plots.}

\item{viol.threshold}{threshold of tolerated constraint violation, used to
determine feasibility if \code{feasible.only == TRUE}.}

\item{nondominated.only}{logical flag to use only nondominated points in the
plots.}

\item{plot.weights}{logical flag to plot the weight vectors for 2 and
3-objective problems.}

\item{which.objectives}{integer vector of which objectives to plot.
Defaults to \code{NULL} (use all objectives)}

\item{suppress.pause}{logical flag to prevent pause messages from being show after every image.
Defaults to \code{FALSE} (show pause messages)}

\item{color.by.obj}{integer, determines which objective is used as the basis
for coloring the parallel coordinates plot.}
}
\description{
S3 method for plotting \emph{moead} objects (the output of \code{\link[=moead]{moead()}}).
}
\section{References}{

F. Campelo, L.S. Batista, C. Aranha (2020): The {MOEADr} Package: A
Component-Based Framework for Multiobjective Evolutionary Algorithms Based on
Decomposition. Journal of Statistical Software \doi{10.18637/jss.v092.i06}\cr
}

\examples{
problem.1 <- list(name = "example_problem",
                  xmin = rep(-1,30),
                  xmax = rep(1,30),
                  m    = 2)
out <- moead(preset    = preset_moead("original2"),
             problem   = problem.1,
             stopcrit  = list(list(name = "maxiter",
                                   maxiter = 100)),
             showpars  = list(show.iters = "dots",
                              showevery  = 10))
plot(out, suppress.pause = TRUE)

}
