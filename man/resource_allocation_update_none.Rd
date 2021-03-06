% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resource_allocation_update_none.R
\name{resource_allocation_update_none}
\alias{resource_allocation_update_none}
\title{Resource Allocation none}
\usage{
resource_allocation_update_none(
  iter,
  resource.allocation,
  priority.values,
  bigZ,
  dt.bigZ,
  neighbors.T,
  Y,
  dt.Y,
  W,
  X,
  dt.X,
  ...
)
}
\arguments{
\item{iter}{integer iteration}

\item{resource.allocation}{list resource allocation method and specific vars}

\item{priority.values}{priority values for deciding which subproblem to give resources}

\item{bigZ}{agg values for RI}

\item{dt.bigZ}{agg values for RI from dt previous iteration}

\item{neighbors.T}{neighborhood size (unused?)}

\item{Y}{PF}

\item{dt.Y}{PF from dt previous iteration}

\item{X}{decision values}

\item{dt.X}{decision values from dt previous iteration}

\item{...}{other parameters (included for compatibility with generic call)}

\item{idx.boundary}{indexes of the boundaries subproblem}

\item{problem}{list with number of objectives n}
}
\value{
priority values for each subproblem \code{out}'.
}
\description{
Proxy Resource Allocation none- dummy
}
\details{
Returns a vector with indexes of each subproblem that will receive resource
}
\section{References}{


F. Campelo, L.S. Batista, C. Aranha (2020): The {MOEADr} Package: A
Component-Based Framework for Multiobjective Evolutionary Algorithms Based on
Decomposition. Journal of Statistical Software \doi{10.18637/jss.v092.i06}\cr
}

