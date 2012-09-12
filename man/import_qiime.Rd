\name{import_qiime}
\alias{import_qiime}
\title{Import function to read files created by the QIIME pipeline.}
\usage{
  import_qiime(otufilename=NULL, mapfilename=NULL,
  treefilename=NULL, biotaxonomy=NULL, showProgress=TRUE,
  chunk.size=1000L, ...)
}
\arguments{
  \item{otufilename}{(Optional). A character string
  indicating the file location of the OTU file. The
  combined OTU abundance and taxonomic identification file,
  tab-delimited, as produced by QIIME under default output
  settings. Default value is \code{NULL}.}

  \item{mapfilename}{(Optional). The QIIME map file
  required for processing pyrosequencing tags in QIIME as
  well as some of the post-clustering analysis. This is a
  required input file for running QIIME. Its strict
  formatting specification should be followed for correct
  parsing by this function. Default value is \code{NULL}.}

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

  \item{biotaxonomy}{(Optional). A character vector
  indicating the name of each taxonomic level in the
  taxonomy-portion of the otu-file, which may not specify
  these levels explicitly. Default value is \code{NULL}.}

  \item{showProgress}{(Optional). A logical. Indicates
  whether import progress/status should be printed to the
  terminal. Default value is \code{TRUE}, meaning the
  progress will be shown.}

  \item{chunk.size}{(Optional). Positive integer. Default
  is \code{1000L}. This is the number of lines to be read
  and processed at a time, passed directly to the
  \code{\link{import_qiime_otu_tax}} function. A lower
  value helps control memory-fault, but slows down the
  import.}

  \item{...}{Additional arguments passed to
  \code{\link{read_tree}}}
}
\value{
  A \code{\link{phyloseq-class}} object.
}
\description{
  QIIME produces several files that can be analyzed in the
  phyloseq-package, including especially an OTU file that
  typically contains both OTU-abundance and taxonomic
  identity information. The map-file is also an important
  input to QIIME that stores sample covariates, converted
  naturally to the \code{\link{sample_data-class}}
  component data type in the phyloseq-package. QIIME may
  also produce a phylogenetic tree with a tip for each OTU,
  which can also be imported by this function.
}
\details{
  See \url{"http://www.qiime.org/"} for details on using
  QIIME. While there are many complex dependencies, QIIME
  can be downloaded as a pre-installed linux virtual
  machine that runs ``off the shelf''.

  The different files useful for import to \emph{phyloseq}
  are not collocated in a typical run of the QIIME
  pipeline. See the main \emph{phyloseq} vignette for an
  example of where ot find the relevant files in the output
  directory.
}
\examples{
# import_qiime(myOtuTaxFilePath, myMapFilePath)
}
\references{
  \url{http://qiime.org/}

  ``QIIME allows analysis of high-throughput community
  sequencing data.'' J Gregory Caporaso, Justin Kuczynski,
  Jesse Stombaugh, Kyle Bittinger, Frederic D Bushman,
  Elizabeth K Costello, Noah Fierer, Antonio Gonzalez Pena,
  Julia K Goodrich, Jeffrey I Gordon, Gavin A Huttley,
  Scott T Kelley, Dan Knights, Jeremy E Koenig, Ruth E Ley,
  Catherine A Lozupone, Daniel McDonald, Brian D Muegge,
  Meg Pirrung, Jens Reeder, Joel R Sevinsky, Peter J
  Turnbaugh, William A Walters, Jeremy Widmann, Tanya
  Yatsunenko, Jesse Zaneveld and Rob Knight; Nature
  Methods, 2010; doi:10.1038/nmeth.f.303
}
\seealso{
  \code{\link{phyloseq}}

  \code{\link{merge_phyloseq}}

  \code{\link{read_tree}}
}

