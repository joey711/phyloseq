################################################################################
#' Convert phyloseq data to MetagenomeSeq MRexperiment object
#'
#' No testing is performed by this function. The phyloseq data is converted
#' to the relevant \code{\link[metagenomeSeq]{MRexperiment-class}} object, which can then be
#' tested in the zero-inflated mixture model framework
#' (e.g. \code{\link[metagenomeSeq]{fitZig}}) 
#' in the metagenomeSeq package.
#' See the
#' \href{http://joey711.github.io/phyloseq-extensions}{phyloseq-extensions}
#' tutorials for more details.
#'
#' @param physeq (Required). \code{\link{phyloseq-class}}.
#' @param ... (Optional). Additional named arguments passed 
#'  to \code{\link[metagenomeSeq]{newMRexperiment}}.
#'  Most users will not need to pass any additional arguments here.
#'  
#' @return A \code{\link[metagenomeSeq]{MRexperiment-class}} object.
#' 
#' @seealso
#'
#'  \code{\link[metagenomeSeq]{fitTimeSeries}}
#'  \code{\link[metagenomeSeq]{fitLogNormal}}
#'  \code{\link[metagenomeSeq]{fitZig}}
#'  \code{\link[metagenomeSeq]{MRtable}}
#'  \code{\link[metagenomeSeq]{MRfulltable}}
#'
#' @export
#' @importFrom Biobase AnnotatedDataFrame
#'  
#' @examples
#'  # Check out the vignette metagenomeSeq for more details.
#'  # vignette("metagenomeSeq")
#'  data(soilrep)
#'  phyloseq_to_metagenomeSeq(soilrep)
phyloseq_to_metagenomeSeq = function(physeq, ...){
  # Enforce orientation. Samples are columns
  if( !taxa_are_rows(physeq) ){ physeq <- t(physeq)}
  # Coerce count data to vanilla matrix of integers
  countData = round(as(otu_table(physeq), "matrix"), digits=0)
  # Create sample annotation if possible
  if(!is.null(sample_data(physeq,FALSE))){
    ADF = AnnotatedDataFrame(data.frame(sample_data(physeq)))  
  } else { 
    ADF = NULL 
  }
  # Create taxa annotation if possible
  if(!is.null(tax_table(physeq,FALSE))){
    TDF = AnnotatedDataFrame(data.frame(OTUname = taxa_names(physeq),
          data.frame(tax_table(physeq)),row.names = taxa_names(physeq)))
  } else {
    TDF = AnnotatedDataFrame(data.frame(OTUname = taxa_names(physeq),
          row.names = taxa_names(physeq)))
  }
  # Create MRexperiment
  if(requireNamespace("metagenomeSeq")){
    mrobj = metagenomeSeq::newMRexperiment(counts = countData, phenoData = ADF, featureData = TDF,...)
    # Calculate normalization factor
    if (sum(colSums(countData > 0) > 1) < ncol(countData)) {
        p = suppressMessages(metagenomeSeq::cumNormStat(mrobj))
    }
    else {
        p = suppressMessages(metagenomeSeq::cumNormStatFast(mrobj))
    }
    mrobj = metagenomeSeq::cumNorm(mrobj, p = p)
    return(mrobj)
  }
}
