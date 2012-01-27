\name{import_mothur_tree}
\alias{import_mothur_tree}
\title{Import and prune mothur-produced tree.}
\usage{
  import_mothur_tree(mothur_tree_file, mothur_list_file,
    cutoff=NULL)
}
\arguments{
  \item{mothur_tree_file}{The tree file name that was
  output from \emph{mothur}.  Probably a file that ends
  with the suffix \code{".tree"}.}

  \item{mothur_list_file}{The list file name and/or
  location as produced by \emph{mothur}.}

  \item{cutoff}{A character string indicating the cutoff
  value, (or \code{"unique"}), that matches one of the
  cutoff-values used to produce the OTU clustering results
  contained within the list-file created by \emph{mothur}
  (and specified by the \code{mothur_list_file} argument).
  The default is to take the largest value among the cutoff
  values contained in the list file. If only one cutoff is
  included in the file, it is taken and this argument does
  not need to be specified. Note that the \code{cluster()}
  function within the \emph{mothur} package will often
  produce a list file with multiple cutoff values, even if
  a specific cutoff is specified. It is suggested that you
  check which cutoff values are available in a given list
  file using the \code{\link{show_mothur_list_cutoffs}}
  function.}
}
\value{
  A reduced tree that has been pruned to contain only
  unique taxa, as defined by a particular OTU-clustering
  cutoff in the \emph{mothur} list file.
}
\description{
  The \code{\link[ape]{read.tree}} function is sufficient
  for importing a \emph{mothur}-produced tree into \code{R}
  in \code{"phylo"} format. This function further requires
  a list file as input, and prunes / renames the tree such
  that it is compatible with the associated abundance data.
  The expected tree output from \emph{mothur} will contain
  a tip for each sequence in the original dataset, even if
  those sequences are highly similar or identical. This
  function returns a reduced tree that has been pruned to
  contain only unique taxa, as defined by a particular
  OTU-clustering cutoff in the \emph{mothur} list file.
}
\details{
  This is a user-available module of a more comprehensive
  function for importing output files from the
  \emph{mothur} package, \code{link{import_mothur}}. The
  \code{import_mothur} function is suggested if the goal is
  to import more than just the tree from \emph{mothur}.
}
\seealso{
  \code{\link{import_mothur}}
}
\keyword{internal}

