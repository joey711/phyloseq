\name{filterfunSample}
\alias{filterfunSample}
\title{A sample-wise generic filter function, analogous to filterfun from the
genefilter package.}
\usage{
  filterfunSample(...)
}
\arguments{
  \item{...}{A comma-separated list of functions.}
}
\value{
  An enclosure (function) that itself will return a logical
  vector, according to the functions provided in the
  argument list, evaluated in order. The output of
  filterfunSample is appropriate for the `flist' argument
  to the genefilterSample method.
}
\description{
  A sample-wise generic filter function, analogous to
  filterfun from the genefilter package.
}
\examples{
#
}
\seealso{
  genefilter genefilterSample
}
\keyword{filter}
\keyword{genefilter}
\keyword{trim}

