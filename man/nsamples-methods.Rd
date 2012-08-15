\docType{methods}
\name{nsamples}
\alias{nsamples}
\alias{nsamples,ANY-method}
\alias{nsamples,otu_table-method}
\alias{nsamples,phyloseq-method}
\alias{nsamples,sample_data-method}
\title{Get the number of samples.}
\usage{
  nsamples(physeq)
}
\arguments{
  \item{physeq}{A \code{\link{phyloseq-class}},
  \code{\link{sample_data}}, or
  \code{\link{otu_table-class}}.}
}
\value{
  An integer indicating the total number of samples.
}
\description{
  Get the number of samples.
}
\examples{
#
# # From "picante" package
# data("phylocom")
# tree <- phylocom$phylo
# OTU1 <- otu_table(phylocom$sample, taxa_are_rows=FALSE)
# nsamples(OTU1)
# physeq1 <- phyloseq(OTU1, tree)
# nsamples(physeq1)
}
\seealso{
  \code{\link{taxa_names}}, \code{\link{sample_names}},
  \code{\link{ntaxa}}
}

