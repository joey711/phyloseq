################################################################################
#' Convert phyloseq data to DESeq2 dds object
#'
#' No testing is performed by this function. The phyloseq data is converted
#' to the relevant \code{\link[DESeq2]{DESeqDataSet}} object, which can then be
#' tested in the negative binomial generalized linear model framework
#' of the \code{\link[DESeq2]{DESeq}} function in DESeq2 package.
#' See the
#' \href{http://joey711.github.io/phyloseq-extensions}{phyloseq-extensions}
#' tutorials for more details.
#'
#' @param physeq (Required). \code{\link{phyloseq-class}}.
#'  Must have a \code{\link{sample_data}} component.
#'
#' @param design (Required). A \code{\link{formula}} which specifies the design of the experiment,
#'  taking the form \code{formula(~ x + y + z)}. That is, a formula with right-hand side only.
#'  By default, the functions in this package and DESeq2
#'  will use the last variable in the formula (e.g. \code{z})
#'  for presenting results (fold changes, etc.) and plotting.
#'  When considering your specification of experimental design, you will want to 
#'  re-order the levels so that the \code{NULL} set is first.
#'  For example, the following line of code would ensure that Enterotype 1 is used as the 
#'  reference sample class in tests by setting it to the first of the factor levels
#'  using the \code{\link{relevel}} function:
#'  
#'  \code{sample_data(entill)$Enterotype <- relevel(sample_data(entill)$Enterotype, "1")}
#'  
#' @param ... (Optional). Additional named arguments passed to \code{\link[DESeq2]{DESeqDataSetFromMatrix}}.
#'  Most users will not need to pass any additional arguments here.
#'  Most testing-related options should be provided in 
#'  a following call to \code{\link[DESeq2]{DESeq}}.
#'  
#' @return A \code{\link[DESeq2]{DESeqDataSet}} object.
#' 
#' @seealso
#' 
#' \code{vignette("phyloseq-mixture-models")}
#' 
#' The 
#' \href{http://joey711.github.io/phyloseq-extensions}{phyloseq-extensions}
#' tutorials.
#' 
#'  \code{\link[DESeq2]{DESeq}}
#'  
#'  \code{\link[DESeq2]{results}}
#'
#' @export
#'  
#' @examples
#'  # Check out the vignette phyloseq-mixture-models for more details.
#'  # vignette("phyloseq-mixture-models")
#'  data(soilrep)
#'  phyloseq_to_deseq2(soilrep, ~warmed)
phyloseq_to_deseq2 = function(physeq, design, ...){
  # Need to add check here for missing sample_data
  if( is.null(sample_data(physeq, FALSE)) ){
    stop("There must be sample_data present, for specifying experimental design. See ?phyloseq_to_deseq2")
  }
  # Enforce orientation. Samples are columns
  if( !taxa_are_rows(physeq) ){ physeq <- t(physeq)}
  # Coerce count data to vanilla matrix of integers
  countData = round(as(otu_table(physeq), "matrix"), digits=0)
  colData = data.frame(sample_data(physeq))
  # Create the DESeq data set, dds.
  if(requireNamespace("DESeq2")){
    dds <- DESeq2::DESeqDataSetFromMatrix(countData, colData, design, ...)
    return(dds)
  }
}
################################################################################
