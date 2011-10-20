\name{topk}
\alias{topk}
\title{The most abundant \code{k} taxa}
\usage{
  topk(k, na.rm=TRUE)
}
\arguments{
  \item{k}{An integer, indicating how many of the most
  abundant taxa should be kept.}

  \item{na.rm}{A logical. Should \code{NA}s be removed.
  Default is \code{TRUE}.}
}
\value{
  Returns a function (enclosure) that will return TRUE for
  each element in the most abundant k values.
}
\description{
  The most abundant \code{k} taxa
}

