\name{import_RDP_otu}
\alias{import_RDP_otu}
\title{Import new RDP OTU-table format}
\usage{
  import_RDP_otu(otufile)
}
\arguments{
  \item{otufile}{(Optional). A character string indicating
  the file location of the OTU file, produced/exported
  according to the instructions above.}
}
\value{
  A \code{\link{otuTable-class}} object.
}
\description{
  Recently updated tools on RDP Pyro site make it easier to
  import Pyrosequencing output into R. The modified tool
  ``Cluster To R Formatter'' can take a cluster file
  (generated from RDP Clustering tools) to create a
  community data matrix file for distance cutoff range you
  are interested in. The resulting output file is a
  tab-delimited file containing the number of sequences for
  each sample for each OTU. The OTU header naming
  convention is \code{"OTU_"} followed by the OTU number in
  the cluster file. It pads ``0''s to make the OTU header
  easy to sort. The OTU numbers are not necessarily in
  order.
}
\examples{
otufile <- "http://cloud.github.com/downloads/joey711/phyloseq/rformat_dist_0.03.txt"
ex_otu  <- import_RDP_otu(otufile)
head(t(ex_otu))
}
\seealso{
  An alternative ``cluster'' file importer for RDP results:
  \code{\link{import_RDP_cluster}}

  The main RDP-pyrosequencing website
  \url{http://pyro.cme.msu.edu/index.jsp}
}

