\name{import_qiime_sampleData}
\alias{import_qiime_sampleData}
\title{Import just \code{sampleData} file from QIIME pipeline.}
\usage{
  import_qiime_sampleData(mapfilename)
}
\arguments{
  \item{mapfilename}{(Required). A character string. The
  name of the QIIME map file required for processing
  pyrosequencing tags in QIIME as well as some of the
  post-clustering analysis. This is a required input file
  for running QIIME. Its strict formatting specification is
  expected by this function, do not attempt to modify it
  manually once it has worked properly in QIIME.}
}
\value{
  A \code{sampleData} object.
}
\description{
  QIIME produces several files that can be analyzed in the
  phyloseq-package, This includes the map-file, which is an
  important \emph{input} to QIIME that can also indicate
  sample covariates. It is converted naturally to the
  sampleData component data type in phyloseq-package, based
  on the R data.frame.
}
\details{
  See \code{\link{import_qiime}} for more information about
  QIIME. It is also the suggested function for importing
  QIIME-produced data files.
}
\seealso{
  \code{\link{import_qiime}}, \code{\link{merge_phyloseq}},
  \code{\link{phyloseq}},
  \code{\link{import_qiime_otu_tax}}
}
\keyword{internal}

