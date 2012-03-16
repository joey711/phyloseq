\name{getVariable}
\alias{getVariable}
\title{Get the values for a particular variable in sampleData}
\usage{
  getVariable(physeq, varName)
}
\arguments{
  \item{physeq}{(Required). \code{\link{sampleData-class}},
  or \code{\link{phyloseq-class}}.}

  \item{varName}{(Required). Character string of the
  variable name in \code{sampleData}. Use
  \code{sample.variables(physeq)} for available variables
  in your object.}
}
\value{
  Data. The clas of the data depends on what the contents
  of sampleData.
}
\description{
  This is a simple accessor function for streamlining
  access to values/vectors/factors/etc contained in the
  sampleData.
}
\examples{
# Load the GlobalPatterns dataset into the workspace environment
data(GlobalPatterns)
# Look at the different values for SampleType
getVariable(GlobalPatterns, "SampleType")
}
\seealso{
  getSpecies species.names sample.names getTaxa
  \code{\link{sample.variables}}
}

