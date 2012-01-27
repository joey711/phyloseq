\name{taxglom.internal}
\alias{taxglom.internal}
\title{taxglom core internal function.}
\usage{
  taxglom.internal(physeq, tax, NArm = TRUE,
    bad_empty = c(NA, "", " ", "\t"))
}
\arguments{
  \item{physeq}{the object on which agglomeration is to
  take place.}

  \item{tax}{for this internal function, tax must be a
  character vector. \code{tax} is a vector in the core
  internal function.}
}
\value{
  the agglomerated object. Class matches argument
  \code{physeq}.
}
\description{
  taxglom.internal makes all the glomming happen, and
  delegates the object-handling issues to
  \code{merge_species()}.
}
\keyword{internal}

