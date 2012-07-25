thisVer <- getRversion()
if(!compareVersion(as.character(thisVer), "2.15.0") >=0)
  stop("R v2.15 is required")
# Define dependencies in installer  
cran.nms <- c("ape", "doParallel", "foreach", "ggplot2", "igraph0", "picante", "vegan", "RJSONIO", "plyr")
bioc.nms <- c("multtest", "genefilter")
##  cran
install.packages(cran.nms)
## BioC
source("http://bioconductor.org/biocLite.R")
biocLite(bioc.nms)
## use devtools to install latest version of phyloseq from GitHub
if (!require('devtools'))
  install.packages('devtools')  # needs Curl for RCurl
library("devtools")
pkgs <- list(joey711 = c("phyloseq"))
for (repo in names(pkgs)) {
  for (pkg in pkgs[[repo]]) install_github(pkg, repo)
}


