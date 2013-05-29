\name{plot_clusgap}
\alias{plot_clusgap}
\title{Create a ggplot summary of gap statistic results}
\usage{
  plot_clusgap(clusgap, title = "Gap Statistic results")
}
\arguments{
  \item{clusgap}{(Required). An object of S3 class
  \code{"clusGap"}, basically a list with components. See
  the \code{\link[cluster]{clusGap}} documentation for more
  details. In most cases this will be the output of
  \code{\link{gapstat_ord}}, or
  \code{\link[cluster]{clusGap}} if you called it
  directly.}

  \item{title}{(Optional). Character string. The main title
  for the graphic. Default is \code{"Gap Statistic
  results"}.}
}
\value{
  A \code{\link[ggplot2]{ggplot}} plot object. The rendered
  graphic should be a plot of the gap statistic score
  versus values for \code{k}, the number of clusters.
}
\description{
  Create a ggplot summary of gap statistic results
}
\examples{
# Load data
data(enterotype)
# ordination
exord = ordinate(enterotype, method = "MDS", distance = "bray")
gs = gapstat_ord(exord, verbose=FALSE, method="Tibs2001SEmax")
print(gs, method="Tibs2001SEmax")
plot_clusgap(gs)
# Non-ordination example
gs.pam.RU <- clusGap(ruspini, FUN = pam1, K.max = 8, B = 60)
gs.pam.RU
plot(gs.pam.RU, main = "Gap statistic for the 'ruspini' data")
mtext("k = 4 is best .. and  k = 5  pretty close")
plot_clusgap(gs.pam.RU)
}
\seealso{
  \code{\link{gapstat_ord}}

  \code{\link[cluster]{clusGap}}

  \code{\link[ggplot2]{ggplot}}
}

