\docType{methods}
\name{index_reorder}
\alias{index_reorder}
\alias{index_reorder,phyloseq-method}
\title{Force index order of phyloseq objects}
\usage{
  index_reorder(ps, index_type)
}
\arguments{
  \item{ps}{(Required). A \code{\link{phyloseq-class}}
  instance.}

  \item{index_type}{(Optional). A character string
  specifying the indices to properly order. Supported
  values are \code{c("both", "taxa", "samples")}. Default
  is \code{"both"}, meaning samples and taxa indices will
  be checked/re-ordered.}
}
\description{
  Force index order of phyloseq objects
}
\examples{
## data("GlobalPatterns")
## GP = index_reorder(GlobalPatterns)
}
\keyword{internal}

