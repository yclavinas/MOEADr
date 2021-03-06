% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preset_moead.R
\name{preset_moead}
\alias{preset_moead}
\title{preset_moead}
\usage{
preset_moead(name = NULL)
}
\arguments{
\item{name}{name of the preset to be generated. Use \code{preset_moead()} to obtain
the list of available options.}
}
\value{
List object containing the preset, to be used as an input to \code{\link[=moead]{moead()}};
or, if \code{name == NULL} (the default), returns a logical flag invisibly.
}
\description{
Generate a preset configuration for moead()].
}
\details{
This function returns a list of configuration presets taken from
the literature to be used with the \code{\link[=moead]{moead()}} function in package \code{MOEADr}.

Use these configurations as a starting point. We strongly
recommend that you play around with the particular configurations
(see example).
}
\section{References}{

F. Campelo, L.S. Batista, C. Aranha (2020): The {MOEADr} Package: A
Component-Based Framework for Multiobjective Evolutionary Algorithms Based on
Decomposition. Journal of Statistical Software \doi{10.18637/jss.v092.i06}\cr
}

\examples{

# Generate list of available presets
preset_moead(name = NULL)

\dontrun{
  library(smoof) # < Install package smoof if needed
  ZDT1 <- make_vectorized_smoof(prob.name  = "ZDT1",
                                dimensions = 30)
                                problem   <- list(name       = "ZDT1",
                                                  xmin       = rep(0, 30),
                                                  xmax       = rep(1, 30),
                                                  m          = 2)

  # Get preset configuration for original MOEA/D
  configuration <- preset_moead("original")

  # Modify whatever you fancy:
  stopcrit <- list(list(name = "maxiter", maxiter = 50))
  showpars <- list(show.iters = "dots", showevery  = 10)
  seed     <- 42

  output <- moead(problem  = problem,
                  preset   = configuration,
                  showpars = showpars,
                  stopcrit = stopcrit,
                  seed     = seed)
}

}
