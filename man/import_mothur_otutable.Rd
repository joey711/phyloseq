\name{import_mothur_otutable}
\alias{import_mothur_otutable}
\title{Import mothur list file and return as list object in R.}
\usage{
  import_mothur_otutable(mothur_list_file,
  mothur_group_file, cutoff=NULL)
}
\arguments{
  \item{mothur_list_file}{The list file name and/or
  location as produced by \emph{mothur}.}

  \item{mothur_group_file}{The name/location of the group
  file produced by \emph{mothur}'s \code{make.group()}
  function. It contains information about the sample source
  of individual sequences, necessary for creating a
  species/taxa abundance table (\code{otuTable}). See
  \code{http://www.mothur.org/wiki/Make.group}}

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
  A list, where each element is a character vector of 1 or
  more sequence identifiers, indicating how each sequence
  from the original data is clustered into OTUs by
  \emph{mothur}. Note that in some cases this is highly
  dependent on the choice for \code{cutoff}.
}
\description{
  This is a user-available module of a more comprehensive
  function for importing OTU clustering/abundance data
  using the \emph{mothur} package. The list object returned
  by this function is not immediately useable by other
  \emph{phyloseq} functions, and must be first parsed in
  conjunction with a separate \emph{mothur} \code{"group"}
  file. This function is made accessible to \emph{phyloseq}
  users for troubleshooting and inspection, but the
  \code{link{import_mothur()}} function is suggested if the
  goal is to import the OTU clustering results from
  \emph{mothur} into \emph{phyloseq}.
}
\seealso{
  \code{\link{show_mothur_list_cutoffs}},
  \code{\link{import_mothur_otulist}}
}

