\name{show_mothur_list_cutoffs}
\alias{show_mothur_list_cutoffs}
\title{Show cutoff values available in a mothur list file}
\usage{
  show_mothur_list_cutoffs(mothur_list_file)
}
\arguments{
  \item{mothur_list_file}{The list file name and/or
  location as produced by \emph{mothur}.}
}
\value{
  A character vector of the different cutoff values
  contained in the file. For a given set of arguments to
  the \code{cluster()} command from within \emph{mothur}, a
  number of OTU-clustering results are returned in the same
  list file. The exact cutoff values used by \emph{mothur}
  can vary depending on the input data. This simple
  function returns the cutoffs that were actually included
  in the \emph{mothur} output. This an important extra step
  prior to importing the OTUs with the
  \code{import_mothur_otulist()} function.
}
\description{
  This is a helper function to report back to the user the
  different cutoff values available in a given \emph{list}
  file created by the OTU clustering and analysis package
  called \emph{mothur}
}
\seealso{
  \code{\link{import_mothur}}
}

