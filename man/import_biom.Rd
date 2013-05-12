\name{import_biom}
\alias{import_biom}
\title{Import phyloseq data from biom-format file}
\usage{
  import_biom(BIOMfilename, treefilename=NULL,
  refseqfilename=NULL, refseqFunction=readDNAStringSet,
  refseqArgs=NULL, parseFunction=parse_taxonomy_default,
  parallel=FALSE, version=1.0, ...)
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

  \item{treefilename}{(Optional). A file representing a
  phylogenetic tree or a \code{\link{phylo}} object. Files
  can be NEXUS or Newick format. See
  \code{\link{read_tree}} for more details. If provided,
  the tree should have the same OTUs/tip-labels as the OTUs
  in the other files. Any taxa or samples missing in one of
  the files is removed from all. For the QIIME pipeline
  this tree is typically a tree of the representative 16S
  rRNA sequences from each OTU cluster, with the number of
  leaves/tips equal to the number of taxa/species/OTUs.
  Default value is \code{NULL}. Note that this argument can
  be a tree object (\code{\link[ape]{phylo}}-class) for
  cases where the tree has been --- or needs to be ---
  imported separately.}

  \item{refseqfilename}{(Optional). Default \code{NULL}.
  The file path of the biological sequence file that
  contains at a minimum a sequence for each OTU in the
  dataset. Alternatively, you may provide an
  already-imported \code{\link[Biostrings]{XStringSet}}
  object that satisfies this condition. In either case, the
  \code{\link{names}} of each OTU need to match exactly the
  \code{\link{taxa_names}} of the other components of your
  data. If this is not the case, for example if the data
  file is a FASTA format but contains additional
  information after the OTU name in each sequence header,
  then some additional parsing is necessary, which you can
  either perform separately before calling this function,
  or describe explicitly in a custom function provided in
  the (next) argument, \code{refseqFunction}. Note that the
  \code{\link[Biostrings]{XStringSet}} class can represent
  any arbitrary sequence, including user-defined
  subclasses, but is most-often used to represent RNA, DNA,
  or amino acid sequences. The only constraint is that this
  special list of sequences has exactly one named element
  for each OTU in the dataset.}

  \item{refseqFunction}{(Optional). Default is
  \code{\link[Biostrings]{readDNAStringSet}}, which expects
  to read a fasta-formatted DNA sequence file. If your
  reference sequences for each OTU are amino acid, RNA, or
  something else, then you will need to specify a different
  function here. This is the function used to read the file
  connection provided as the the previous argument,
  \code{refseqfilename}. This argument is ignored if
  \code{refseqfilename} is already a
  \code{\link[Biostrings]{XStringSet}} class.}

  \item{refseqArgs}{(Optional). Default \code{NULL}.
  Additional arguments to \code{refseqFunction}. See
  \code{\link[Biostrings]{XStringSet-io}} for details about
  additional arguments to the standard read functions in
  the Biostrings package.}

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
  \href{http://www.qiime.org/svn_documentation/documentation/biom_format.html}{biom-format}
}
\seealso{
  \code{\link{import}}

  \code{\link{import_qiime}}

  \code{\link{read_tree}}

  \code{\link[biom]{read_biom}}

  \code{\link[biom]{biom_data}}

  \code{\link[biom]{sample_metadata}}

  \code{\link[biom]{observation_metadata}}

  \code{\link[Biostrings]{XStringSet-io}}
}

