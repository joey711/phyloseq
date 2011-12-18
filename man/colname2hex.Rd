\name{colname2hex}
\alias{colname2hex}
\title{Convert color name to hex value.}
\usage{
  colname2hex(colname, alpha = 1)
}
\arguments{
  \item{colname}{A character vector of \code{R} color
  names.}

  \item{alpha}{The standard alpha parameter specifying
  opacity/transparency.}
}
\value{
  A color hex value of length equal to length of
  \code{colname}.
}
\description{
  For converting color name to color hex and providing an
  optional opacity parameter (alpha)
}
\seealso{
  col2rgb rgb2hsv colors
}
\keyword{internal}

