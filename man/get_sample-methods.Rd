\docType{methods}
\name{get_sample}
\alias{get_sample}
\alias{get_sample,otu_table-method}
\alias{get_sample,phyloseq-method}
\title{Returns all abundance values for species \code{i}.}
\usage{
  get_sample(physeq, i)
}
\arguments{
  \item{physeq}{(Required). \code{\link{otu_table-class}},
  or \code{\link{phyloseq-class}}.}

  \item{i}{(Required). A single taxa/species/OTU ID for
  which you want to know the abundance in each sample.}
}
\value{
  An integer vector of the abundance values for each sample
  in \code{physeq} for species \code{i}
}
\description{
  This is a simple accessor function for investigating a
  single species-of-interest.
}
\examples{
data(esophagus)
taxa_names(esophagus)
get_sample(esophagus, "59_5_19")
}
\seealso{
  \code{\link{get_taxa}} \code{\link{taxa_names}}
  \code{\link{sample_names}}
}

