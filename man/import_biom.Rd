\name{import_biom}
\alias{import_biom}
\title{Import phyloseq data from biom-format file}
\usage{
  import_biom(BIOMfilename, tree=NULL,
  parseFunction=parse_taxonomy_default, parallel=FALSE,
  version=0.9, ...)
}
\arguments{
  \item{BIOMfilename}{(Required). A character string
  indicating the file location of the BIOM formatted file.
  This is a JSON formatted file, specific to biological
  datasets, as described in
  \url{http://www.qiime.org/svn_documentation/documentation/biom_format.html}{the
  biom-format home page}. In principle, this file should
  include you OTU abundance data (OTU table), your
  taxonomic classification data (taxonomy table), as well
  as your sample data, for instance what might be in your
  ``sample map'' in QIIME. A phylogenetic tree is not yet
  supported by biom-format, and so is a separate argument
  here. If, for some reason, your biom-format file is
  missing one of these mentioned data types but you have it
  in a separate file, you can first import the data that is
  in the biom file using this function, \code{import_biom},
  and then ``merge'' the remaining data after you have
  imported with other tools using the relatively
  general-purpose data merging function called
  \code{\link{merge_phyloseq}}.}

  \item{tree}{(Optional). Character string or ``phylo''
  class. A file path to the phylogenetic tree file
  associated with the data in the \code{BIOMfilename}, or
  an already-imported version of that phylogenetic tree as
  a phylo object. Import of the tree will be attempted
  using \code{\link{read_tree}}, and additional import
  parameters for that function can be supplied (see
  below).}

  \item{parseFunction}{(Optional). A function. It must be a
  function that takes as its first argument a character
  vector of taxonomic rank labels for a single OTU and
  parses and names each element (an optionally removes
  unwanted elements). Further details and examples of
  acceptable functions are provided in the documentation
  for \code{\link{parse_taxonomy_default}}. There are many
  variations on taxonomic nomenclature, and naming
  conventions used to store that information in various
  taxonomic databases and phylogenetic assignment
  algorithms. A popular database,
  \url{http://greengenes.lbl.gov/cgi-bin/nph-index.cgi}{greengenes},
  has its own custom parsing function provided in the
  phyloseq package,
  \code{\link{parse_taxonomy_greengenes}}, and more can be
  contributed or posted as code snippets as needed. They
  can be custom-defined by a user immediately prior to the
  the call to \code{\link{import_biom}}, and this is a
  suggested first step to take when trouble-shooting
  taxonomy-related errors during file import.}

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

  \item{...}{Additional parameters passed on to
  \code{\link{read_tree}}.}
}
\value{
  A \code{\link{phyloseq-class}} object.
}
\description{
  New versions of QIIME produce a more-comprehensive and
  formally-defined JSON file format, called biom file
  format:
}
\details{
  ``The biom file format (canonically pronounced `biome')
  is designed to be a general-use format for representing
  counts of observations in one or more biological samples.
  BIOM is a recognized standard for the Earth Microbiome
  Project and is a Genomics Standards Consortium candidate
  project.''

  \url{http://biom-format.org/}
}
\examples{
# An included example of a rich dense biom file
rich_dense_biom  <- system.file("extdata", "rich_dense_otu_table.biom",  package="phyloseq")
import_biom(rich_dense_biom,  parseFunction=parse_taxonomy_greengenes)
# An included example of a sparse dense biom file
rich_sparse_biom <- system.file("extdata", "rich_sparse_otu_table.biom", package="phyloseq")
import_biom(rich_sparse_biom, parseFunction=parse_taxonomy_greengenes)
# # # Example code for importing large file with parallel backend
# library("doParallel")
# registerDoParallel(cores=6)
# import_biom("my/file/path/file.biom", parseFunction=parse_taxonomy_greengenes, parallel=TRUE)
}
\references{
  \url{http://www.qiime.org/svn_documentation/documentation/biom_format.html}{biom-format}
}
\seealso{
  \code{\link{import}}

  \code{\link{import_qiime}}

  \code{\link{read_tree}}
}

