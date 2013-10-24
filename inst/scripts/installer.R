thisVer <- getRversion()
if(!compareVersion(as.character(thisVer), "3.0.0") >=0)
  stop("R v3.0 or greater is required")
## Step 1
## Install the release version from BioC
source("http://bioconductor.org/biocLite.R")
biocLite("phyloseq")
## You may also want to go further and install the latest devel from BioC,
## or further still, from GitHub (the most recent).
## Step 2 
## Install from BioC-devel
devel = "http://bioconductor.org/packages/2.14/bioc"
biocLite("phyloseq", siteRepos=devel, suppressUpdates=TRUE, type="source")
## Step 3
## Use devtools::install_github to install latest version of phyloseq from GitHub
if(!require("devtools", quietly=TRUE)){
  # Note: needs Curl for RCurl
  install.packages("devtools")
}
library("devtools")
devtools::install_github("phyloseq", "joey711")
