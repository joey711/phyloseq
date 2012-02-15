\name{import_biom}
\alias{import_biom}
\title{Import phyloseq data from BIOM file}
\usage{
  import_biom(BIOMfilename, taxaPrefix=NULL,
  parallel=FALSE, version=0.9)
}
\arguments{
  \item{BIOMfilename}{(Required). A character string
  indicating the file location of the BIOM formatted file.
  This is a JSON formatted file, specific to biological
  datasets, as described in

  \url{http://www.qiime.org/svn_documentation/documentation/biom_format.html}}

  \item{taxaPrefix}{(Optional). Character string. What
  category of prefix precedes the taxonomic label at each
  taxonomic rank. Currently only ``greengenes'' is a
  supported option, and implies that the first letter
  indicates the taxonomic rank, followed by two underscores
  and then the actual taxonomic assignment at that rank.
  The default value is \code{NULL}, meaning that no prefix
  or rank identifier will be interpreted.}

  \item{parallel}{(Optional). Logical. Wrapper option for
  \code{.parallel} parameter in \code{plyr-package}
  functions. If \code{TRUE}, apply parsing functions in
  parallel, using parallel backend provided by
  \code{\link{foreach}} and its supporting backend
  packages. One caveat, plyr-parallelization currently
  works most-cleanly with \code{multicore}-like backends
  (Mac OS X, Unix?), and may throw warnings for SNOW-like
  backends. See the example below for code invoking
  multicore-style backend within the \code{doParallel}
  package.

  Finally, for many datasets a parallel import should not
  be necessary because a serial import will be just as fast
  and the import is often only performed one time; after
  which the data should be saved as an RData file using the
  \code{\link{save}} function.}

  \item{version}{(Optional). Numeric. The expected version
  number of the file. As the BIOM format evolves,
  version-specific importers will be available by adjusting
  the version value. Default is \code{0.9}. Not
  implemented. Has no effect (yet).}
}
\value{
  A \code{\link{phyloseq-class}} object.
}
\description{
  New versions of QIIME produce a more-comprehensive and
  formally-defined JSON file format. From the QIIME
  website:
}
\details{
  ``The biom file format (canonically pronounced `biome')
  is designed to be a general-use format for representing
  counts of observations in one or more biological
  samples.''

  \url{http://www.qiime.org/svn_documentation/documentation/biom_format.html}
}
\examples{
# # # import with default parameters, specify a file
# import_BIOM(myBIOMfile)
# # # Example code for importing large file with parallel backend
# library("doParallel")
# registerDoParallel(cores=6)
# import_biom("my/file/path/file.biom", taxaPrefix="greengenes", parallel=TRUE)
}
\references{
  \url{http://www.qiime.org/svn_documentation/documentation/biom_format.html}
}
\seealso{
  \code{\link{import}}, \code{\link{import_qiime}}
}

