\docType{methods}
\name{import}
\alias{import}
\alias{import,character-method}
\title{Universal import method (wrapper) for phyloseq-package}
\usage{
  import(pipelineName, ...)
}
\arguments{
  \item{pipelineName}{(Required). Character string. The
  name of the analysis tool / pipeline / package that
  created the OTU-cluster data or other data that you now
  want to import. Current options are \code{c("mothur",
  "pyrotagger", "QIIME", "RDP")}, and only the first letter
  is necessary.}

  \item{...}{(Required). Additional arguments providing
  file paths, and possible other paramaters to the desired
  tool-specific import function.}
}
\value{
  In most cases a \code{\link{phyloseq-class}} will be
  returned, though the type of individual component data
  types will vary according to the type(s) of data created
  by the respective tool, and whether those files were
  provided as arguments. The expected behavior of the
  import methods is to return as comprehensive an object as
  possible, containing all the available data types of the
  arguments provided.
}
\description{
  A user must still understand the additional arguments
  required for each type of import data. Those arguments
  are described in detail at the tool-specific
  \code{import_*} links below. Each clustering tool /
  package / pipeline has its own idiosyncratic set of file
  names / types, and it remains the responsibility of the
  user to understand which file-path should be provided to
  each argument for the particular importing submethod.
  This method merely provides a central documentation and
  method-name, and the arguments are passed along as-is.
}
\examples{
## import("QIIME", otufilename=myOtuTaxFilePath, mapfilename=myMapFilePath)
}
\references{
  mothur: \url{http://www.mothur.org/wiki/Main_Page}

  PyroTagger: \url{http://pyrotagger.jgi-psf.org/}

  QIIME: \url{http://qiime.org/}

  RDP pipeline: \url{http://pyro.cme.msu.edu/index.jsp}
}
\seealso{
  For mothur, see: \code{\link{import_mothur}}

  Separate tools for mothur are also:
  \code{\link{show_mothur_list_cutoffs}}
  \code{\link{import_mothur_dist}}
  \code{\link{export_mothur_dist}}

  For PyroTagger, see: \code{\link{import_pyrotagger_tab}}

  For QIIME, see: \code{\link{import_qiime}}

  For RDP pipeline, see: \code{\link{import_RDP_cluster}}
}

