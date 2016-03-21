################################################################################
#' Convert phyloseq OTU count data into DGEList for edgeR package
#' 
#' @param physeq (Required).  A \code{\link{phyloseq-class}} or
#'  an \code{\link{otu_table-class}} object. 
#'  The latter is only appropriate if \code{group} argument is also a 
#'  vector or factor with length equal to \code{nsamples(physeq)}.
#'  
#' @param group (Required). A character vector or factor giving the experimental
#'  group/condition for each sample/library. Alternatively, you may provide
#'  the name of a sample variable. This name should be among the output of
#'  \code{sample_variables(physeq)}, in which case
#'  \code{get_variable(physeq, group)} would return either a character vector or factor.
#'  This is passed on to \code{\link[edgeR]{DGEList}}, 
#'  and you may find further details or examples in its documentation.
#'  
#' @param method (Optional). The label of the edgeR-implemented normalization to use.
#'  See \code{\link[edgeR]{calcNormFactors}} for supported options and details. 
#'  The default option is \code{'RLE'}, which is a scaling factor method 
#'  proposed by Anders and Huber (2010).
#'  At time of writing, the \link[edgeR]{edgeR} package supported 
#'  the following options to the \code{method} argument:
#'  
#'  \code{c('TMM', 'RLE', 'upperquartile', 'none')}.
#'
#' @param ... Additional arguments passed on to \code{\link[edgeR]{DGEList}}
#' 
#' @export
#' 
#' @importFrom edgeR DGEList
#' @importFrom edgeR estimateTagwiseDisp
#' @importFrom edgeR estimateCommonDisp
#' @importFrom edgeR calcNormFactors

phyloseq_to_edgeR <- function(physeq, group, method = "RLE", ...) {
  # Enforce orientation.
  if (!taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }
  x <- as(otu_table(physeq), "matrix")
  # Add one to protect against overflow, log(0) issues.
  x <- x + 1
  # Check `group` argument
  if (identical(all.equal(length(group), 1), TRUE) & nsamples(physeq) > 1) {
    # Assume that group was a sample variable name (must be categorical)
    group <- get_variable(physeq, group)
  }
  # Define gene annotations (`genes`) as tax_table
  taxonomy <- tax_table(physeq, errorIfNULL = FALSE)
  if (!is.null(taxonomy)) {
    taxonomy <- data.frame(as(taxonomy, "matrix"))
  }
  # Now turn into a DGEList
  y <- DGEList(counts = x, group = group, genes = taxonomy, remove.zeros = TRUE, ...)
  # Calculate the normalization factors
  z <- calcNormFactors(y, method = method)
  # Check for division by zero inside `calcNormFactors`
  if (!all(is.finite(z$samples$norm.factors))) {
    stop("Something wrong with edgeR::calcNormFactors on this data, \n 
         non-finite $norm.factors, consider changing `method` argument")
  }
  # Estimate dispersions
  return(estimateTagwiseDisp(estimateCommonDisp(z)))
  }


################################################################################
#' Plot edgeR results for a phyloseq or a edgeR object.
#' @param data (Required): a \code{\link{phyloseq-class}} or a \code{\link{DESeqDataSet-class}} object.
#'#' @param tax_table: Require if data is a \code{\link{DESeqDataSet-class}} object. 
#' The taxonomic table used to find the \code{taxa} and \code{color_taxa} arguments. If data is a \code{\link{phyloseq-class}} object, data@tax_table is used.
#' @param contrast (Required):This argument specifies what comparison to extract from the object to build a results table. See \code{\link[DESeq2]{results}} man page for more details.
#' @param alpha (default = 0.01): the significance cutoff used for optimizing the independent filtering. If the adjusted p-value cutoff (FDR) will be a value other than 0.1, alpha should be set to that value.
#' @param taxa (default = 'Genus'): taxonomic level of interest
#' @param color_tax (default = 'Phylum'): taxonomic level used for color assignation
#' @param verbose : whether the function print some information during the computation
#' @param ... Additional arguments passed on to \code{\link[edgeR]{exactTest}} or \code{\link[ggplot2]{ggplot}}
#' 
#' @export
#' 
#' @examples 
#'\dontrun{
#'data(GlobalPatterns)
#'plot_edgeR_phyloseq(GlobalPatterns, c('SampleType', 'Soil', 'Feces'), tax_table = GlobalPatterns@tax_table, color_tax = 'Kingdom')
#'plot_edgeR_phyloseq(GlobalPatterns, c('SampleType', 'Soil', 'Feces'), taxa = 'Class', tax_table = GlobalPatterns@tax_table, color_tax = 'Kingdom')
#'}
#' @author Adrien Taudiere
#'
#' @return A \code{\link{ggplot}}2 plot representing edgeR results
#'
#' @seealso \code{\link[edgeR]{exactTest}}
#' @seealso \code{\link{plot_deseq2_phyloseq}}

plot_edgeR_phyloseq <- function(data, contrast = NULL, alpha = 0.01, taxa = "Genus", color_tax = "Phylum", 
                                verbose = TRUE, ...) {
  
  if (!inherits(data, "phyloseq")) {
    stop("data must be an object of class 'phyloseq'")
  }
  
  if (verbose) {
    message("Conversion to edgeR format")
  }
  data_edgeR <- phyloseq_to_edgeR(data, group = contrast[1])
  if (verbose) {
    message("Perform edgeR binary test")
  }
  et <- exactTest(data_edgeR, pair = c(contrast[2], contrast[3]), ...)
  
  tt <- topTags(et, n = nrow(et$table), adjust.method = "BH", sort.by = "PValue")
  res <- tt@.Data[[1]]
  sigtab <- res[(res$FDR < alpha), ]
  sigtab <- cbind(as(sigtab, "data.frame"))
  
  sigtabgen <- subset(sigtab, !is.na(taxa))
  
  d <- tapply(sigtabgen$logFC, sigtabgen[, color_tax], function(x) max(x))
  d <- sort(d, TRUE)
  sigtabgen$col_tax <- factor(as.character(sigtabgen[, color_tax]), levels = names(d))
  
  d <- tapply(sigtabgen$logFC, sigtabgen[, taxa], function(x) max(x))
  d <- sort(d, TRUE)
  sigtabgen$tax <- factor(as.character(sigtabgen[, taxa]), levels = names(d))
  
  p <- ggplot(sigtabgen, aes(x = tax, y = logFC, color = col_tax), ...) + geom_point(size = 6) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) + 
    labs(title = paste("Change in abundance for ", contrast[1], " (", contrast[2], " vs ", contrast[3], ")", sep = ""))
  
  print(p)
  
  return(invisible(p))
}
################################################################################

