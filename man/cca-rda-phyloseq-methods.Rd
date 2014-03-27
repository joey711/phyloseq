\docType{methods}
\name{cca.phyloseq}
\alias{cca.phyloseq}
\alias{cca.phyloseq,otu_table-method}
\alias{cca.phyloseq,phyloseq,formula-method}
\alias{cca.phyloseq,phyloseq,NULL-method}
\alias{rda.phyloseq}
\title{Constrained Correspondence Analysis and Redundancy Analysis.}
\arguments{
  \item{physeq}{(Required). Phylogenetic sequencing data
  (\code{\link{phyloseq-class}}). The data on which you
  want to perform the ordination.}

  \item{formula}{(Optional). A \code{\link{formula}},
  specifying the contraining variable(s) format, with
  variable names corresponding to \code{\link{sample_data}}
  (RHS) from within \code{physeq}.}

  \item{method}{(Optional). A single
  \code{\link{character}} string, specifying \code{"RDA"}
  or \code{"CCA"}. Default is \code{"CCA"}.}

  \item{...}{(Optional). Additional named arguments passed
  to \code{\link[vegan]{capscale}}.}
}
\value{
  same output as \code{\link[vegan]{cca}} or
  \code{\link[vegan]{rda}}, respectively.
}
\description{
  This is the internal function that simplifies getting
  phyloseq data into the constrained ordination functions,
  \code{\link[vegan]{cca}} and \code{\link[vegan]{rda}}.
  Unlike \code{\link[phyloseq]{capscale.phyloseq}}, the
  formula argument to these methods is optional, and
  results in an unconstrained ordination.
}
\examples{
#
# cca.phyloseq(physeq, formula, method, ...)
}
\seealso{
  \code{\link{plot_ordination}}, \code{\link[vegan]{rda}},
  \code{\link[vegan]{cca}}
}
\keyword{internal}

