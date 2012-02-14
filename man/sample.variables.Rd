\name{sample.variables}
\alias{sample.variables}
\title{Get the sample variables present in sampleData}
\usage{
  sample.variables(physeq)
}
\arguments{
  \item{physeq}{(Required). \code{\link{sampleData-class}},
  or \code{\link{phyloseq-class}}.}
}
\value{
  Character vector. The names of the variables in the
  sampleData data.frame. Essentially the column names.
  Useful for selecting model and graphics parameters that
  interact with sampleData.
}
\description{
  This is a simple accessor function to make it more
  convenient to determine the sample variable names of a
  particular \code{\link{phyloseq-class}} object.
}
\examples{
data(enterotype)
sample.variables(enterotype)
}
\seealso{
  getSpecies species.names sample.names getTaxa
}

