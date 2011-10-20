\name{import_qiime_sampleMap}
\alias{import_qiime_sampleMap}
\title{Import just OTU/Taxonomy file from QIIME pipeline.}
\usage{
  import_qiime_sampleMap(mapfilename)
}
\arguments{
  \item{mapfilename}{The QIIME map file required for
  processing pyrosequencing tags in QIIME as well as some
  of the post-clustering analysis. This is a required input
  file for running QIIME. Its strict formatting
  specification should be followed for correct parsing by
  this function.}
}
\value{
  A \code{sampleMap} object.
}
\description{
  QIIME produces several files that can be analyzed in the
  phyloseq-package, This includes the map-file, which is an
  important \emph{input} to QIIME that can also indicate
  sample covariates. It is converted naturally to the
  sampleMap component data type in phyloseq-package, based
  on the R data.frame.
}
\details{
  Add reference to the QIIME pipeline.
}
\examples{
#
}
\seealso{
  \code{\link{import_qiime}}, \code{\link{merge_phyloseq}},
  \code{\link{phyloseq}},
  \code{\link{import_qiime_otu_tax}}
}

