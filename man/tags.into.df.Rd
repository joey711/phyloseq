\name{tags.into.df}
\alias{tags.into.df}
\title{Convert tree tags into data.frame columns}
\usage{
  tags.into.df(phylo, df)
}
\arguments{
  \item{phylo,}{input phylo object}

  \item{df,}{data.frame with a 'node' column corresponding
  to integer indices of tree nodes.}
}
\value{
  df, a copy of the input data frame, with tags added as
  new columns
}
\description{
  Given a \code{\link{phylo}} object and a data frame,
  transform all the tags from the tree into columns of the
  data frame. Rows of the data frame are linked to the tree
  via a required 'node' column, which must contain integer
  indices of the associated node.
}
\details{
  This function is similar to the tree.as.data.frame
  method, but not exactly the same. It is used internally
  by the tree.layout function.
}
\author{
  Gregory Jordan \email{gjuggler@gmail.com}
}
\keyword{internal}

