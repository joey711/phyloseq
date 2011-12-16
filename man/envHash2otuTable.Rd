\name{envHash2otuTable}
\alias{envHash2otuTable}
\title{Convert a sequence-sample hash (like ENV file) into an OTU table.}
\usage{
  envHash2otuTable(tipSampleTable)
}
\arguments{
  \item{tipSampleTable}{(Required). A two-column character
  table (matrix or data.frame), where each row specifies
  the sequence name and source sample, consistent with the
  env-file for the UniFrac server
  (\url{http://bmf2.colorado.edu/unifrac/}).}
}
\value{
  \code{\link{otuTable}}. A trivial OTU table where each
  sequence is treated as a separate OTU.
}
\description{
  Parses an ENV-file into a sparse matrix of
  species-by-sample, where each species-row has only one
  non-zero value. We call this sparse abundance table the
  trivial OTU table, where every sequence is treated as a
  separate species. If a phylogenetic tree is available, it
  can be submitted with this table as arguments to
  \code{\link{tipglom}} to create an object with a
  non-trivial \code{otuTable}.
}
\examples{
#
## fakeSeqNameVec <- paste("seq_", 1:8, sep="")
## fakeSamNameVec <- c(rep("A", 4), rep("B", 4))
## fakeSeqAbunVec <- sample(1:50, 8, TRUE)
## test    <- cbind(fakeSeqNameVec, fakeSamNameVec, fakeSeqAbunVec)
## testotu <- envHash2otuTable( test )
## test    <- cbind(fakeSeqNameVec, fakeSamNameVec)
## testotu <- envHash2otuTable( test )
}
\references{
  \url{http://bmf2.colorado.edu/unifrac/}
}
\seealso{
  \code{\link{import_env_file}}, \code{\link{tipglom}},
  \code{\link{otuTable}}
}
\keyword{internal}

