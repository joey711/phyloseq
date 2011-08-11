\name{read_env_file}
\alias{read_env_file}
\title{Read a UniFrac-formatted ENV file.}
\usage{
  read_env_file(envfilename, tree = NULL, sep = "\t")
}
\arguments{
  \item{envfilename}{A charater string of the ENV filename
  (relative or absolute)}

  \item{tree}{Optional phylo object that will be used to
  prune elements in the ENV hash table once it is read.
  Taxa not present as tip.labels in \code{tree} will be
  excluded from the output.}

  \item{sep}{A character string indicating the delimiter
  used in the file. The default is \code{"\t"}.}
}
\value{
  A two- or three- column hash table (matrix).
}
\description{
  Convenience wrapper function to read the
  environment-file, as formatted for input to the UniFrac
  server(REF). The official format of these files is that
  each row specifies (in order) the sequence name, source
  sample, and (optionally) the number of times the sequence
  was observed.
}
\examples{
#
}
\keyword{environment}
\keyword{file}
\keyword{OTU}
\keyword{trivial}

