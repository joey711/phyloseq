\name{import_qiime_otu_tax}
\alias{import_qiime_otu_tax}
\title{Import just OTU/Taxonomy file from QIIME pipeline.}
\usage{
  import_qiime_otu_tax(otufilename, biotaxonomy = NULL)
}
\arguments{
  \item{otufilename}{A character string indicating the file
  location of the OTU file. The combined OTU abundance and
  taxonomic identification file, tab-delimited, as produced
  by QIIME under default output settings.  Default value is
  \code{NULL}.}

  \item{biotaxonomy}{A character vector indicating the name
  of each taxonomic level in the taxonomy-portion of the
  otu-file, which may not specify these levels explicitly.
  Default value is \code{NULL}.}
}
\value{
  An \code{otuTax} object.
}
\description{
  QIIME produces several files that can be analyzed in the
  phyloseq-package, including especially an OTU file that
  typically contains both OTU-abundance and taxonomic
  identity information.
}
\details{
  Add reference to the QIIME pipeline.
}
\seealso{
  \code{\link{merge_phyloseq}}, \code{\link{phyloseq}},
  \code{\link{import_qiime_sampleData}}
}
\keyword{internal}

