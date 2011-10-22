\name{parapply}
\alias{parapply}
\title{Wrapper for parallelization}
\usage{
  parapply(i, fun1, num_cores, ...)
}
\arguments{
  \item{i}{(Required). The vector/list to loop through.}

  \item{fun1}{(Required). The function to apply.}

  \item{num_cores}{(Required). Integer. The number of cores
  to use.}

  \item{...}{(Optional). Additional arguments, passed to
  fun1.}
}
\value{
  a sample-by-sample distance matrix, suitable for NMDS,
  etc.
}
\description{
  This function is a convenience wrapper around
  \code{Rmpi::mpi.spawn.Rslaves()} and especially
  \code{Rmpi::mpi.apply}. In essence, it should act as a
  parallel apply-loop on the iteration vector, \code{i}.
}
\keyword{internal}

