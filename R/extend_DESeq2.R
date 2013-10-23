################################################################################
#' Convert phyloseq data to DESeq2 dds object
#'
#' No testing is performed by this function. The phyloseq data is converted
#' to the relevant \code{\link[DESeq2]{DESeqDataSet}} object, which can then be
#' tested in the negative binomial generalized linear model framework
#' of the \code{\link[DESeq2]{DESeq}} function in DESeq2 package.
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
#'  \code{\link[DESeq2]{DESeq}}
#'  
#'  \code{\link[DESeq2]{results}}
#'
#' @export
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#'  
#' @examples
#'  ### Load and process data.
#' data("enterotype")
#' # Transform back to counts from samples, 
#' # assume that the smallest non-zero value, snzv, in each sample is equal to 1.
#' entcount = transform_sample_counts(enterotype, function(x){
#'    snzv = min(x[x > 0])
#'    return(round(x/snzv))
#' })
#' # Remove the odd "-1" and "Bacteria" garbage genera that includes all the unannotated sequences
#' # (almost certainly a mixture of many genera and noise)
#' entcount = prune_taxa(!taxa_names(entcount) %in% c("-1", "Bacteria"), entcount)
#' # Subset to just illumina samples
#' entill = subset_samples(entcount, SeqTech == "Illumina")
#' dds = phyloseq_to_deseq2(entill, ~ Enterotype)
#' # The DESeq() and results() methods come from DESeq2-package
#' library("DESeq2")
#' dds = DESeq(dds, fitType="local")
#' # The default multiple-inference correction is BH, already occurred in DESeq
#' res = results(dds)
#' res = res[order(res$padj, na.last=NA), ]
#' # Print the results for OTUs with adjusted-P above 0.05
#' print(res[(res$padj < 0.05), ])
phyloseq_to_deseq2 = function(physeq, design, ...){
  # Need to add check here for missing sample_data
  if( is.null(sample_data(physeq, FALSE)) ){
    stop("There must be sample_data present, for specifying experimental design. See ?phyloseq_to_deseq2")
  }
  # Enforce orientation. Samples are columns
  if( !taxa_are_rows(physeq) ){ physeq <- t(physeq)}
  # Coerce count data to vanilla matrix of integers
  countData = round(as(otu_table(physeq), "matrix"), digits=0)
  countData = countData + 1L
  colData = data.frame(sample_data(physeq))
  # Create the DESeq data set, dds.
  dds <- DESeqDataSetFromMatrix(countData, colData, design, ...)
  return(dds)
}
################################################################################
