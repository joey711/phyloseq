##############################
# An example function for installing phyloseq from various sources
##############################
install_phyloseq = function(branch = "release",
                            minRVersion = "3.3.0",
                            verbose = TRUE){
  if(!compareVersion(as.character(getRversion()), minRVersion) >=0){
    stop("phyloseq installation script failed.\n", 
         "R ", minRVersion, " or greater is required.")
  }
  branch <- as.character(branch)
  if(branch == "release"){
    if(verbose){
      message("Installing the release version from BioC")
    }
    source("http://bioconductor.org/biocLite.R")
    biocLite("phyloseq", suppressUpdates=TRUE)
    return("release")
  }
  if(branch == "devel"){
    if(verbose){
      message("Installing the devel version from BioC")
    }
    biocLite("phyloseq", 
             siteRepos="http://bioconductor.org/packages/devel/bioc",
             suppressUpdates=TRUE,
             type="source")
    return("devel")
  }
  if(branch == "github"){
    if(verbose){
      message("Installing the devel version from joey711/mater from GitHub")
    }
    if(!require("devtools", quietly=TRUE)){
      # Note: needs Curl for RCurl
      install.packages("devtools")
    }
    library("devtools")
    devtools::install_github("phyloseq", "joey711")
    return("github")
  }
  return("something didn't work well")
}
###############
# Execute the function w/ default params.
# You can select alternatives if you want :-)
###############
install_phyloseq()
