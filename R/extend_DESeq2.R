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
################################################################################
################################################################################
################################################################################



####################################################################################
#Plot the result of a DESeq2 test
####################################################################################
#' Plot DESeq2 results for a phyloseq or a DESeq2 object.
#' @param data (Required): a \code{\link{phyloseq-class}} or a \code{\link{DESeqDataSet-class}} object.
#' @param tax_table : Required if data is a \code{\link{DESeqDataSet-class}} object. 
#' The taxonomic table used to find the \code{taxa} and \code{color_taxa} arguments. If data is a \code{\link{phyloseq-class}} object, data@tax_table is used.
#' @param contrast (Required):This argument specifies what comparison to extract from the object to build a results table. See \code{\link[DESeq2]{results}} man page for more details.
#' @param alpha (default = 0.01): the significance cutoff used for optimizing the independent filtering. If the adjusted p-value cutoff (FDR) will be a value other than 0.1, alpha should be set to that value.
#' @param taxa (default = 'Genus'): taxonomic level of interest
#' @param color_tax (default = 'Phylum'): taxonomic level used for color or a color vector.
#' @param taxDepth (default = NULL): Taxonomic depth to test for differential distribution among contrast. If Null the analysis is done at the OTU (i.e. Species) level. If not Null data need to be a \code{\link{phyloseq-class}} object.
#' @param verbose : whether the function print some information during the computation
#' @param ... Additional arguments passed on to \code{\link[DESeq2]{DESeq}} or \code{\link[ggplot2]{ggplot}}
#'
#' @importFrom DESeq2 DESeq
#' @importFrom DESeq2 results
#' @export
#'
#' @examples
#'\dontrun{
#'data("GlobalPatterns")
#'GlobalPatterns_Proteobacteria <- subset_taxa(GlobalPatterns, GlobalPatterns@tax_table[, 2] == 'Proteobacteria')
#'res_deseq2 <- DESeq(phyloseq_to_deseq2(GlobalPatterns_Proteobacteria, ~ SampleType), test = 'Wald', fitType = 'local')
#'plot_deseq2_phyloseq(res_deseq2, c('SampleType', 'Soil', 'Feces'), tax_table = GlobalPatterns_Proteobacteria@tax_table, color_tax = 'Kingdom')
#'plot_deseq2_phyloseq(res_deseq2, c('SampleType', 'Soil', 'Feces'), tax_table = GlobalPatterns_Proteobacteria@tax_table, color_tax = c("red", "black"), alpha = 0.7))
#'plot_deseq2_phyloseq(GlobalPatterns_Proteobacteria, c('SampleType', 'Soil', 'Feces'), taxDepth = 'Family', taxa = 'Family',  color_tax = 'Class')
#'}
#' @author Adrien Taudiere
#' 
#' @return A \code{\link{ggplot}}2 plot representing DESeq2 results
#'
#' @seealso \code{\link[DESeq2]{DESeq}}
#' @seealso \code{\link[DESeq2]{results}}
#' @seealso \code{\link{plot_edgeR_phyloseq}}

plot_deseq2_phyloseq <- function(data, contrast = NULL, tax_table = NULL, alpha = 0.01, 
                                 taxa = "Genus", color_tax = "Phylum", taxDepth = NULL, verbose = TRUE, ...) {
  
  if (!inherits(data, "phyloseq")) {
    if (!inherits(data, "DESeqDataSet")) {
      stop("data must be an object of class 'phyloseq' or 'DESeqDataSet'")
    }
  } else {
    # Calculate new dataset given the Taxa depth if taxDepth is not null
    if (!is.null(taxDepth)) {
      data_TAX <- data
      data_TAX@otu_table <- otu_table(apply(data@otu_table, 2, function(x) tapply(x, 
                                                                                  data@tax_table[, taxDepth], sum)), taxa_are_rows = T)
      data_TAX@tax_table <- tax_table(apply(data@tax_table[, 1:match(taxDepth, colnames(data@tax_table))], 
                                            2, function(x) X <- tapply(x, data@tax_table[, taxDepth], function(xx) xx[1])))
      data_TAX@refseq <- NULL
      data <- data_TAX
      if (is.na(match(taxa, colnames(data@tax_table)))) {
        taxa <- taxDepth
      }
    }
    
    if (is.null(tax_table) & inherits(data, "phyloseq")) {
      tax_table <- data@tax_table
    }
    
    if (verbose) {
      message("Conversion to Deseq2 format.")
    }
    data_deseq2 <- phyloseq_to_deseq2(data, as.formula(paste("~", contrast[1])))
    
    if (verbose) {
      message("Calculation of Deseq2 results.")
    }
    data <- DESeq(data_deseq2, test = "Wald", fitType = "parametric", quiet = !verbose, 
                  ...)
  }
  
  # Calcul deseq2 results
  res <- results(data, contrast = contrast)
  
  d <- res[which(res$padj < alpha), ]
  
  if (dim(d)[1] == 0) {
    message("None taxa present significant distribution pattern through contrast.")
    return("None taxa present significant distribution pattern through contrast.")
  }
  d <- cbind(as(d, "data.frame"), as(tax_table[rownames(d), ], "matrix"))
  
  # Compute colors
  areColors <- function(x) {
    sapply(x, function(X) {
      tryCatch(is.matrix(col2rgb(X)), error = function(e) FALSE)
    })
  }
  
  if (!sum(areColors(color_tax)) > 0) {
    x <- tapply(d$log2FoldChange, d[, color_tax], function(x) max(x))
    x <- sort(x, TRUE)
    d$col_tax <- factor(as.character(d[, color_tax]), levels = names(x))
  } else {
    d$col_tax <- rep(color_tax, length = dim(d)[1])
  }
  
  # Compute log2FoldChange values
  x <- tapply(d$log2FoldChange, d[, taxa], function(x) max(x))
  x <- sort(x, TRUE)
  d$tax <- factor(as.character(d[, taxa]), levels = names(x))
  
  if (!sum(areColors(color_tax)) > 0) {
    p <- ggplot(d, aes(x = tax, y = log2FoldChange, color = col_tax), ...) + geom_point(size = 6) + 
      theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) + labs(title = paste("Change in abundance for ", 
                                                                                                  contrast[1], " (", contrast[2], " vs ", contrast[3], ")", sep = ""))
  } else {
    p <- ggplot(d, aes(x = tax, y = log2FoldChange), ...) + geom_point(size = 6, color = d$col_tax) + 
      theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) + labs(title = paste("Change in abundance for ", 
                                                                                                  contrast[1], " (", contrast[2], " vs ", contrast[3], ")", sep = ""))
  }
  print(p)
  
  return(invisible(p))
}

