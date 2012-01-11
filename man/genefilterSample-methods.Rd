\docType{methods}
\name{genefilterSample}
\alias{genefilterSample}
\alias{genefilterSample,matrix-method}
\alias{genefilterSample,otuTable-method}
\alias{genefilterSample,phyloseq-method}
\title{Filter OTUs with arbitrary function, sample-wise.}
\arguments{
  \item{X}{The object that needs trimming. Can be matrix,
  otuTable, or higher- order phyloseq classes that contain
  an otuTable.}

  \item{flist}{An enclosure object, typically created with
  \code{\link{filterfunSample}}}

  \item{A}{An integer. The number of samples in which a
  taxa / species passed the filter for it to be labeled
  TRUE in the output logical vector.}
}
\value{
  A logical vector with names equal to species.names (or
  rownames, if matrix).
}
\description{
  A general OTU trimming function for selecting OTUs that
  satisfy some criteria within the distribution of each
  sample, and then also an additional criteria for number
  of samples that must pass. This is a genefilter-like
  function that only considers sample-wise criteria. The
  number of acceptable samples is used as the final
  criteria (set by the argument \code{A}) to determine
  whether or not the taxa should be retained (\code{TRUE})
  or not (\code{FALSE}). Just like with genefilter, a
  logical having length equal to
  nrow()/\code{\link{nspecies}} is returned, indicating
  which should be kept. This output can be provided
  directly to OTU trimming function,
  \code{\link{prune_species}}. By contrast,
  \code{\link[genefilter]{genefilter}}, of the genefilter
  package in Bioconductor, works only on the rows of a
  matrix. Note that, because \code{\link{otuTable-class}}
  inherits directly from the \code{\link{matrix-class}}, an
  unmodified otuTable can be provided to \code{genefilter},
  but be mindful of the orientation of the otuTable (use
  \code{\link{speciesAreRows}}), and transpose
  (\code{\link[phyloseq]{t}}) if needed.
}
\examples{
#
## testOTU <- otuTable(matrix(sample(1:50, 25, replace=TRUE), 5, 5), speciesAreRows=FALSE)
## f1  <- filterfunSample(topk(2))
## wh1 <- genefilterSample(testOTU, f1, A=2)
## wh2 <- c(T, T, T, F, F)
## prune_species(wh1, testOTU)
## prune_species(wh2, testOTU)
##
## taxtab1 <- taxTab(matrix("abc", 5, 5))
## prune_species(wh1, taxtab1)
## prune_species(wh2, taxtab1)
}
\seealso{
  \code{\link[genefilter]{genefilter}},
  \code{\link{filterfunSample}}, \code{\link[phyloseq]{t}},
  \code{\link{prune_species}}
}
\keyword{agglomerate}
\keyword{cluster}
\keyword{OTU}
\keyword{tree}

