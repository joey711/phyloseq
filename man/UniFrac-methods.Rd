\docType{methods}
\name{UniFrac}
\alias{UniFrac}
\alias{UniFrac,phyloseq-method}
\title{Calculate weighted or unweighted UniFrac distance for all sample pairs.}
\arguments{
  \item{physeq}{(Required). \code{\link{phyloseq-class}},
  containing at minimum a phylogenetic tree
  (\code{\link{phylo-class}}) and contingency table
  (\code{\link{otuTable-class}}). See examples below for
  coercions that might be necessary.}

  \item{weighted}{(Optional). Logical. Should use
  weighted-UniFrac calculation?  Weighted-UniFrac takes
  into account the relative abundance of species/taxa
  shared between samples, whereas unweighted-UniFrac only
  considers presence/absence. Default is \code{FALSE},
  meaning the unweighted-UniFrac distance is calculated for
  all pairs of samples.}

  \item{normalized}{(Optional). Logical. Should the output
  be normalized such that values range from 0 to 1
  independent of branch length values? Default is
  \code{TRUE}.  Note that (unweighted) \code{UniFrac} is
  always normalized by total branch-length, and so this
  value is ignored when \code{weighted == FALSE}.}

  \item{parallel}{(Optional). Logical. Should execute
  calculation in parallel, using multiple CPU cores
  simultaneously? This can dramatically hasten the
  computation time for this function. However, it also
  requires that the user has registered a parallel
  ``backend'' prior to calling this function.  Default is
  \code{FALSE}. If FALSE, UniFrac will register a serial
  backend so that \code{foreach::\%dopar\%} does not throw
  a warning.}
}
\value{
  a sample-by-sample distance matrix, suitable for NMDS,
  etc.
}
\description{
  This function calculates the UniFrac distance for all
  sample-pairs in a \code{\link{phyloseq-class}} object
  using the abundances (\code{\link{otuTable-class}}) and a
  phylogenetic tree (\code{\link{phylo-class}}). If the
  tree and contingency table are separate objects,
  suggested solution is to combine them using the
  \code{\link{phyloseq}} function. For example, the
  following code
}
\details{
  \code{phyloseq(myOTUtable, myTree)}

  returns a \code{phyloseq}-class object that has been
  pruned and comprises the minimum arguments necessary for
  \code{UniFrac()}.

  Parallelization is possible, and encouraged for this
  computing-intensive calculation. It has been implemented
  with the \emph{foreach} package. This means that parallel
  calls need to be preceded by 2 or more commands that
  register the parallel ``backend''. This is acheived via
  your choice of helper packages. One of the simplest seems
  to be the \emph{doMC} package. See the commented-out
  examples at the bottom for running \code{UniFrac()} in
  parallel via doMC. Note that doMC was chosen as an
  example because it is relatively easy to use, but it does
  not work for Windows (or other non-posix) systems. For
  Windows, try ``doSNOW'' or ``doParallel'' as an adaptor
  package for parallelized UniFrac calculations.

  For more information, see the following links on
  registering the ``backend'':

  \emph{foreach} package manual:

  \url{http://cran.r-project.org/web/packages/foreach/index.html}

  Notes on parallel computing in \code{R}. Skip to the
  section describing the \emph{foreach Framework}. It gives
  off-the-shelf examples for registering a parallel backend
  using the \emph{doMC}, \emph{doSNOW}, or \emph{doMPI}
  packages:

  \url{http://trg.apbionet.org/euasiagrid/docs/parallelR.notes.pdf}

  Furthermore, as of \code{R} version \code{2.14.0} and
  higher, a parallel package is included as part of the
  core installation, \code{\link{parallel-package}}, and
  this can be used as the parallel backend with the
  \code{\link{foreach-package}} using the adaptor package
  ``doParallel''.

  \url{http://cran.r-project.org/web/packages/doParallel/index.html}
}
\examples{
# ################################################################################
# # Perform UniFrac on esophagus data
# ################################################################################
# data("esophagus")
# (y <- UniFrac(esophagus, TRUE))
# UniFrac(esophagus, TRUE, FALSE)
# UniFrac(esophagus, FALSE)
# picante::unifrac(as(t(otuTable(esophagus)), "matrix"), tre(esophagus))
# ################################################################################
# # Try phylocom example data from picante package
# # It comes as a list, so you must construct the phyloseq object first.
# ################################################################################
# data("phylocom")
# otu  <- otuTable(phylocom$sample, FALSE)
# (x1   <- phyloseq(otu, phylocom$phylo))
# UniFrac(x1, TRUE)
# UniFrac(x1, TRUE, FALSE)
# UniFrac(x1, FALSE)
# picante::unifrac(phylocom$sample, phylocom$phylo)
##
}
\references{
  \url{http://bmf.colorado.edu/unifrac/}

  Lozupone, Hamady and Knight, ``UniFrac - An Online Tool
  for Comparing Microbial Community Diversity in a
  Phylogenetic Context.'', BMC Bioinformatics 2006, 7:371

  Lozupone, Hamady, Kelley and Knight, ``Quantitative and
  qualitative (beta) diversity measures lead to different
  insights into factors that structure microbial
  communities.'' Appl Environ Microbiol. 2007

  Lozupone C, Knight R. ``UniFrac: a new phylogenetic
  method for comparing microbial communities.'' Appl
  Environ Microbiol. 2005 71 (12):8228-35.
}
\seealso{
  \code{vegan::vegdist}
}

