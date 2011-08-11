\name{rda.phyloseq}
\alias{rda.phyloseq}
\title{phyloseq package extension for \code{vegan::rda}.}
\arguments{
  \item{X}{formula-class object specifying the input. No
  need to splat higher-order objects. \code{rda.phyloseq}
  understands where to find the abundance table and sample
  characteristics.}

  \item{data}{\code{data.frame} containing information
  equivalent to a \code{sampleMap} object / component. Only
  necessary if complex object does not already contain
  \code{sampleMap} or you are keeping the data separate for
  some reason.}
}
\value{
  same output as \code{\link{rda}}
}
\description{
  A formula is main input to \code{vegan::cca}. This
  complicates dispatch based on object signature. A new
  method with a separate name is defined instead.
}
\examples{
#
## data(ex1)
## rda.phyloseq(ex1 ~ Diet + Gender)
}
\seealso{
  cca.phyloseq
}
\keyword{analysis}
\keyword{rda}
\keyword{RDA}
\keyword{redundancy}

