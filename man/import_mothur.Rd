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
\examples{
# # The following example assumes you have downloaded the esophagus example
# # dataset from the mothur wiki:
# # "http://www.mothur.org/wiki/Esophageal_community_analysis"
# # "http://www.mothur.org/w/images/5/55/Esophagus.zip"
# # The path on your machine may (probably will) vary
# mothur_list_file  <- "~/Downloads/mothur/Esophagus/esophagus.an.list"
# mothur_group_file <- "~/Downloads/mothur/Esophagus/esophagus.good.groups"
# mothur_tree_file  <- "~/Downloads/mothur/Esophagus/esophagus.tree"
# # # Actual examples follow:
# show_mothur_list_cutoffs(mothur_list_file)
# test1 <- import_mothur(mothur_list_file, mothur_group_file, mothur_tree_file)
# test2 <- import_mothur(mothur_list_file, mothur_group_file, mothur_tree_file, cutoff="0.02")
# # Returns just a tree
# import_mothur(mothur_list_file, mothur_tree_file=mothur_tree_file)
# # Returns just an otuTable
# import_mothur(mothur_list_file, mothur_group_file=mothur_group_file)
# # Returns an error
# import_mothur(mothur_list_file)
# # Should return an "OMG, you must provide the list file" error
# import_mothur()
}
\references{
  \url{http://www.mothur.org/wiki/Main_Page}

  Schloss, P.D., et al., Introducing mothur: Open-source,
  platform-independent, community-supported software for
  describing and comparing microbial communities. Appl
  Environ Microbiol, 2009. 75(23):7537-41.
}

