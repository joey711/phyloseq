\name{import_mothur}
\alias{import_mothur}
\title{General function for importing mothur files into phyloseq.}
\usage{
  import_mothur(mothur_list_file, mothur_group_file=NULL,
  mothur_tree_file=NULL, cutoff=NULL)
}
\arguments{
  \item{mothur_list_file}{Required. The list file name /
  location produced by \emph{mothur}.}

  \item{mothur_group_file}{Optional. The name/location of
  the group file produced by \emph{mothur}'s
  \code{make.group()} function. It contains information
  about the sample source of individual sequences,
  necessary for creating a species/taxa abundance table
  (\code{otuTable}). See
  \code{http://www.mothur.org/wiki/Make.group}}

  \item{mothur_tree_file}{Optional. The tree file name
  produced by \emph{mothur}. Probably a file that ends with
  the suffix \code{".tree"}.}

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
  The object class depends on the provided arguments. If
  the first three arguments are provided, then an
  \code{otuTree} object should be returned, containing both
  an OTU-only tree and its associated \code{otuTable}-class
  abundance table. If only a list and group file are
  provided, then an \code{otuTable} object is returned.
  Similarly, if only a list and tree file are provided,
  then only a tree is returned (\code{"phylo"} class).
}
\description{
  General function for importing mothur files into
  phyloseq.
}

