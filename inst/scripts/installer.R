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
    if (!requireNamespace("BiocManager", quietly=TRUE))
        install.packages("BiocManager")
    BiocManager::install("phyloseq", suppressUpdates=TRUE)
    return("phyloseq installed from BioC release branch (if no errors).")
  }
  if(branch == "devel"){
    if(verbose){
      message("\n\nInstalling phyloseq from the devel version from BioC...\n")
    }
    BiocManager::install("phyloseq", 
             siteRepos="http://bioconductor.org/packages/devel/bioc",
             suppressUpdates=TRUE,
             type="source")
    return("phyloseq installed from BioC devel branch (if no errors).")
  }
  if(branch == "github"){
    if(verbose){
      message("Installing the devel version from joey711/master from GitHub")
    }
    if(!require("devtools", quietly=TRUE)){
      # Note: needs Curl for RCurl
      install.packages("devtools")
    }
    library("devtools")
    devtools::install_github("joey711/phyloseq")
    return("phyloseq installed from GitHub `joey711/phyloseq` (if no errors).")
  }
  return("You probably selected an unsupported argument to `branch`.
  Try again using 'release', 'devel', or 'github'.")
}
###############
# Execute the function w/ default params.
# You can select alternatives if you want :-)
###############
install_phyloseq()
