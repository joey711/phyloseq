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

  \item{parallel}{(Optional). Logical. Wrapper option over
  \code{.parallel} parameter in \code{plyr-package}
  functions. If TRUE, apply parsing functions in parallel,
  using parallel backend provided by \code{\link{foreach}}
  and its supporting backen packages.}

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
# import_BIOM(myBIOMfile)
}
\references{
  \url{http://www.qiime.org/svn_documentation/documentation/biom_format.html}
}
\seealso{
  \code{\link{import}}, \code{\link{import_qiime}}
}

