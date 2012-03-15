\name{reconcile_categories}
\alias{reconcile_categories}
\title{Cleans absent levels in sampleData/data.frame.}
\usage{
  reconcile_categories(DFSM)
}
\arguments{
  \item{DFSM}{(Required). A \code{data.frame} or
  \code{sampleData} object that needs to be cleaned.}
}
\value{
  A single \code{data.frame} object. Even if the input
  argument is a \code{sampleData}, the return is a
  \code{data.frame}. Because this is intended to be used
  internally by the builder method, it cannot also call the
  builder function to re-build the cleaned
  \code{sampleData}.
}
\description{
  This is used internally by the builder method,
  \code{\link{sampleData}}, to ensure that the factors
  describing categorical variables in a data.frame or
  sampleData object are free of extra levels that can
  plague downstream plots analysis.
}
\examples{
# # # data(GlobalPatterns)
# # # SM <- sampleData(GlobalPatterns)
# # # DF <- data.frame(SM)
# # # DF <- data.frame(DF, col1=1:nrow(DF), col2=paste(1:nrow(DF), "t", sep=""))
# # # DF <- reconcile_variables(DF)
# # # SM <- sampleData(reconcile_variables(SM))
# # # sapply(DF, class)
# # # sapply(SM, class)
}
\keyword{internal}

