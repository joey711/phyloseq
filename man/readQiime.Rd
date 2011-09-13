\name{readQiime}
\alias{readQiime}
\title{Import function to read output from the QIIME pipeline.}
\usage{
  readQiime(otufilename = NULL, mapfilename = NULL,
  treefilename = NULL, biotaxonomy = NULL, ...)
}
\arguments{
  \item{otufilename}{A character string indicating the file
  location of the OTU file. The combined OTU abundance and
  taxonomic identification file, tab-delimited, as produced
  by QIIME under default output settings. Default value is
  \code{NULL}.}

  \item{mapfilename}{The QIIME map file required for
  processing pyrosequencing tags in QIIME as well as some
  of the post-clustering analysis. This is a required input
  file for running QIIME. Its strict formatting
  specification should be followed for correct parsing by
  this function. Default value is \code{NULL}.}

  \item{treefilename}{A phylogenetic tree in NEXUS format.
  For the QIIME pipeline this is typically a tree of the
  representative 16S rRNA sequences from each OTU cluster,
  with the number of leaves/tips equal to the number of
  taxa/species/OTUs. Default value is \code{NULL}.
  ALTERNATIVELY, this argument can be a tree object
  ()\code{"phylo4"} or \code{"phylo"} class), in case the
  tree has already been imported, or is in a different
  format than NEXUS.}

  \item{biotaxonomy}{A character vector indicating the name
  of each taxonomic level in the taxonomy-portion of the
  otu-file, which may not specify these levels explicitly.
  Default value is \code{NULL}.}

  \item{...}{Additional arguments passed to
  \code{\link{read.nexus}}, as necessary. Make sure that
  your phylogenetic tree file is readable by
  \code{\link{read.nexus}} prior to calling this function.}
}
\value{
  The class of the object returned by \code{readQiime}
  depends upon which filenames are provided. The most
  comprehensive class is chosen automatically, based on the
  input files listed as arguments. \code{readQiime()} will
  return nothing.
}
\description{
  QIIME produces several files that can be analyzed in the
  phyloseq-package, including especially an OTU file that
  typically contains both OTU-abundance and taxonomic
  identity information. The map-file is also an important
  input to QIIME that can also indicate sample covariates,
  converted naturally to the sampleMap component data type
  in phyloseq-package.
}
\details{
  Add reference to the QIIME pipeline.
}
\examples{
#
}
\seealso{
  \code{\link{merge_phyloseq}}, \code{\link{phyloseq}},
  \code{\link{readQiime_otu_tax}},
  \code{\link{readQiime_sampleMap}},
  \code{\link{read.tree}}, \code{\link{read.nexus}},
  \code{\link{readNexus}} \code{\link{readNewick}}
}

