\docType{methods}
\name{merge_samples}
\alias{merge_samples}
\alias{merge_samples,otuTable-method}
\alias{merge_samples,phyloseq-method}
\alias{merge_samples,sampleData-method}
\title{Merge samples based on a sample variable or factor.}
\usage{
  merge_samples(x, group, fun=mean)
}
\arguments{
  \item{x}{(Required). An instance of a phyloseq class that
  has sample indices. This includes
  \code{\link{sampleData-class}},
  \code{\link{otuTable-class}}, and
  \code{\link{phyloseq-class}}.}

  \item{group}{(Required). Either the a single character
  string matching a variable name in the corresponding
  sampleData of \code{x}, or a factor with the same length
  as the number of samples in \code{x}.}

  \item{fun}{(Optional). The function that will be used to
  merge the values that correspond to the same group for
  each variable. It must take a numeric vector as first
  argument and return a single value. Default is
  \code{\link[base]{mean}}. Note that this is (currently)
  ignored for the otuTable, where the equivalent function
  is \code{\link[base]{sum}}, but evaluated via
  \code{\link[base]{rowsum}} for efficiency.}
}
\value{
  A phyloseq object that has had its sample indices merged
  according to the factor indicated by the \code{group}
  argument. The output class matches \code{x}.
}
\description{
  The purpose of this method is to merge/agglomerate the
  sample indices of a phyloseq object according to a
  categorical variable contained in a sampleData or a
  provided factor.
}
\details{
  NOTE: (\code{\link[ape]{phylo}}) trees and
  \code{\link{taxonomyTable-class}} are not modified by
  this function, but returned in the output object as-is.
}
\examples{
#
# data(GlobalPatterns)
# t1 <- merge_samples(sampleData(GlobalPatterns), "SampleType")
# t4 <- merge_samples(GlobalPatterns, "SampleType")
# identical(t1, sampleData(t4))
}
\seealso{
  \code{\link{merge_species}}, code{\link{merge_phyloseq}}
}

