% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resource_allocation_init.R
\name{resource_allocation_init}
\alias{resource_allocation_init}
\title{Resource Allocation Init}
\usage{
resource_allocation_init(resource.allocation, W, bigZ, ...)
}
\arguments{
\item{resource.allocation}{which Resource Allocation method}

\item{W}{matrix of weight}

\item{bigZ}{bigZ legado}

\item{...}{other parameters (included for compatibility with generic call)}
}
\value{
list of data  \code{out}'.
}
\description{
Resource Allocation Initialiation
}
\details{
Create initial structures for Resource Allocation in MOEA/D-ps
indexes and its boundaries
}
\section{References}{


F. Campelo, L.S. Batista, C. Aranha (2020): The {MOEADr} Package: A
Component-Based Framework for Multiobjective Evolutionary Algorithms Based on
Decomposition. Journal of Statistical Software \doi{10.18637/jss.v092.i06}\cr
}

