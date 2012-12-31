\name{import_qiime_otu_tax}
\alias{import_qiime_otu_tax}
\title{Import a QIIME-formatted otu-tax file into a list of two matrices.}
\usage{
  import_qiime_otu_tax(file,
  parseFunction=parse_taxonomy_qiime, parallel=FALSE,
  chunk.size=1000, verbose=TRUE)
}
\arguments{
  \item{file}{(Required). The path to the qiime-formatted
  file you want to import into R. Can be compressed (e.g.
  \code{.gz}, etc.), though the details may be OS-specific.
  That is, Windows-beware.}

  \item{parseFunction}{(Optional). An optional custom
  function for parsing the character string that contains
  the taxonomic assignment of each OTU. The default parsing
  function is \code{\link{parse_taxonomy_qiime}},
  specialized for splitting the \code{";"}-delimited
  strings and also attempting to interpret greengenes
  prefixes, if any, as that is a common format of the
  taxonomy string produced by QIIME.}

  \item{parallel}{(Optional). Logical. Should the parsing
  be performed in parallel?. Default is \code{FALSE}. Only
  a few steps are actually parallelized, and for most
  datasets it will actually be faster and more efficient to
  keep this set to \code{FALSE}. Also, to get any benefit
  at all, you will need to register a parallel ``backend''
  through one of the backend packages supported by the
  \code{\link{foreach-package}}.}

  \item{chunk.size}{(Optional). Positive integer. Default
  is \code{1000L}. This is the number of lines to be read
  and processed at a time. In many cases this will have a
  much larger effect on speed and efficiency than the
  setting for \code{parallel}. If you set this value too
  large, you risk a memory fault. If you set it too low,
  you risk a very slow parsing process.}

  \item{verbose}{(Optional). Logical. Default is
  \code{TRUE}. Should status updates of the parsing process
  be printed to screen?}
}
\value{
  A list of two matrices. \code{$otutab} contains the OTU
  Table as a numeric matrix, while \code{$tax_table}
  contains a character matrix of the taxonomy assignments.
}
\description{
  QIIME produces several files that can be analyzed in the
  phyloseq-package, including especially an OTU file that
  typically contains both OTU-abundance and taxonomic
  identity information.
}
\details{
  See \url{"http://www.qiime.org/"} for details on using
  QIIME. While there are many complex dependencies, QIIME
  can be downloaded as a pre-installed linux virtual
  machine that runs ``off the shelf''.

  This function uses chunking to perform both the reading
  and parsing in blocks of optional size, thus constrain
  the peak memory usage. feature should make this importer
  accessible to machines with modest memory, but with the
  caveat that the full numeric matrix must be a manageable
  size at the end, too. In principle, the final tables will
  be large, but much more efficiently represented than the
  character-stored numbers. If total memory for storing the
  numeric matrix becomes problematic, a switch to a sparse
  matrix representation of the abundance -- which is
  typically well-suited to this data -- might provide a
  solution.
}
\examples{
otufile <- system.file("extdata", "GP_otu_table_rand_short.txt.gz", package="phyloseq")
import_qiime_otu_tax(otufile)
}
\seealso{
  \code{\link{import}}, \code{\link{merge_phyloseq}},
  \code{\link{phyloseq}}, \code{\link{import_qiime}}
  \code{\link{import_qiime_otu_tax}}
  \code{\link{import_env_file}}
}

