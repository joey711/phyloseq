\docType{methods}
\name{UniFrac}
\alias{UniFrac}
\alias{UniFrac,otuTable,phylo-method}
\alias{UniFrac,otuTree,ANY-method}
\title{Calculate weighted or unweighted UniFrac for all samples in an OTU table.}
\usage{
  UniFrac(OTU, tree, weighted=FALSE, normalized=TRUE,
  parallel=FALSE)
}
\arguments{
  \item{OTU}{(Required). \code{otuTable}, or an
  \code{otuTree}. If you have instead a simple matrix of
  abundances, see \code{\link{otuTable}} for coercing it to
  the \code{otuTable} class.}

  \item{tree}{(Optional). Object of class \code{phylo}. Not
  necessary (and ignored) if \code{OTU} is an object that
  also contains a tree (\code{otuTree} class, or its
  children).}

  \item{weighted}{(Optional). Logical. Should use
  weighted-UniFrac calculation? Weighted-UniFrac takes into
  account the relative abundance of species/taxa shared
  between samples, whereas unweighted-UniFrac only
  considers presence/absence. Default is \code{FALSE},
  meaning the unweighted-UniFrac distance is calculated for
  all pairs of samples.}

  \item{normalized}{(Optional). Logical. Should the output
  be normalized such that values range from 0 to 1
  independent of branch length values? Default is
  \code{TRUE}. Note that (unweighted) \code{UniFrac} is
  always normalized by total branch-length, and so this
  value is ignored when \code{weighted == FALSE}.}

  \item{parallel}{(Optional). Logical. Should execute
  calculation in parallel, using multiple CPU cores
  simultaneously? This can dramatically hasten the
  computation time for this function. However, it also
  requires that the user has registered a parallel
  ``backend'' prior to calling this function. Default is
  \code{FALSE}. If FALSE, UniFrac will register a serial
  backend so that \code{foreach::\%dopar\%} does not throw
  a warning.}
}
\value{
  a sample-by-sample distance matrix, suitable for NMDS,
  etc.
}
\description{
  This function calculates the weighted-UniFrac distance
  for all sample-pairs in a species-abundance table using
  the abundances and a phylogenetic tree. If \code{OTU} is
  a more complex object that already contains a
  phylogenetic tree and abundance table, then the argument
  \code{tree} is not necessary and will be ignored.
  WARNING: The species names of both \code{OTU} and
  \code{tree} must match exactly, or weird index erros can
  result. A warning has been added to further protect users
  from encountering this issue unknowingly. The easiest
  solution is to combine the \code{OTU} and \code{tree}
  objects using the \code{\link{phyloseq}} function. E.g.
  \code{phyloseq(OTU, tree)} will return an
  \code{otuTree}-class object that has been pruned and
  comprises the minimum arguments necessary for
  \code{UniFrac()}.
}
\details{
  Parallelization is possible, and encouraged for this
  computing-intensive calculation. It has been implemented
  with the \emph{foreach} package. This means that parallel
  calls need to be preceded by 2 or more commands that that
  register the parallel ``backend''. This is acheived via
  your choice of helper packages. One of the simplest seems
  to be the \emph{doMC} package. See the commented-out
  examples at the bottom for running \code{UniFrac()} in
  parallel via doMC.

  For more information, see the following links on
  registering the ``backend'':

  \emph{foreach} package manual:

  \code{http://cran.r-project.org/web/packages/foreach/index.html}

  Notes on parallel computing in \code{R}. Skip to the
  section describing the \emph{foreach Framework}. It gives
  off-the-shelf examples for registering a parallel backend
  using the \emph{doMC}, \emph{doSNOW}, or \emph{doMPI}
  packages:

  \code{http://trg.apbionet.org/euasiagrid/docs/parallelR.notes.pdf}
}
\examples{
#
}
\seealso{
  \code{vegan::vegdist}
}

