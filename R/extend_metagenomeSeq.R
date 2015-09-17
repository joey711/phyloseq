################################################################################
#' Convert phyloseq data to MetagenomeSeq MRexperiment object
#'
#' No testing is performed by this function. The phyloseq data is converted
#' to the relevant \code{\link[metagenomeSeq]{'MRexperiment-class'}} object, which can then be
#' tested in the zero-inflated mixture model framework
#' of the \code{\link[metagenomeSeq]{fitFeatureModel}} or \code{\link[metagenomeSeq]{fitZig}} functions in the metagenomeSeq package.
#' See the
#' \href{http://joey711.github.io/phyloseq-extensions}{phyloseq-extensions}
#' tutorials for more details.
#'
#' @param physeq (Required). \code{\link{phyloseq-class}}.
#'  
#' @param ... (Optional). Additional named arguments passed to \code{\link[metagenomeSeq]{newMRexperiment}}.
#'  Most users will not need to pass any additional arguments here.
#'  Most testing-related options should be provided in 
#'  a following call to \code{\link[metagenomeSeq]{fitFeatureModel}}.
#'  
#' @return A \code{\link[metagenomeSeq]{'MRexperiment-class'}} object.
#' 
#' @seealso
#'
#'  \code{\link[metagenomeSeq]{fitTimeSeries}}
#'
#'  \code{\link[metagenomeSeq]{fitLogNormal}}
#'
#'  \code{\link[metagenomeSeq]{fitFeatureModel}}
#'
#'  \code{\link[metagenomeSeq]{fitZig}}
#' 
#'  \code{\link[metagenomeSeq]{MRtable}}
#'  
#'  \code{\link[metagenomeSeq]{MRfulltable}}
#'
#' @export
#'
#' @importFrom metagenomeSeq newMRexperiment
#'
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
  OTU = as(otu_table(physeq), "matrix")
  if(!is.null(sample_data(physeq,FALSE))){
    ADF = AnnotatedDataFrame(data.frame(sample_data(physeq)))  
  }
  if(!is.null(tax_table(physeq,FALSE))){
    TDF = AnnotatedDataFrame(data.frame(OTUname = taxa_names(physeq),
          data.frame(tax_table(physeq)),row.names = taxa_names(physeq)))
  }
  mrobj = newMRexperiment(counts = OTU, phenoData = ADF, featureData = TDF)
  if (sum(colSums(OTU > 0) > 1) < ncol(OTU)) {
      p = suppressMessages(cumNormStat(mrobj))
  }
  else {
      p = suppressMessages(cumNormStatFast(mrobj))
  }
  mrobj = cumNorm(mrobj, p = p)
  return(mrobj)
}
