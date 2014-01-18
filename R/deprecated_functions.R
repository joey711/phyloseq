################################################################################
#' Depcrecated functions in the phyloseq package.
#' 
#' These will be migrated to \code{"defunct"} status in the next release,
#' and removed completely in the release after that.
#' These functions are provided for compatibility with older version of
#' the phyloseq package.  They may eventually be completely
#' removed.
#' 
#' @usage deprecated_phyloseq_function(x, value, ...)
#' @rdname phyloseq-deprecated
#' @name phyloseq-deprecated
#' @param x For assignment operators, the object that will undergo a replacement
#'  (object inside parenthesis).
#' @param value For assignment operators, the value to replace with 
#'  (the right side of the assignment).
#' @param ... For functions other than assignment operators, 
#'  parameters to be passed to the modern version of the function (see table).
#' @docType package
#' @export plot_taxa_bar taxaplot taxtab taxTab sampleData samData sam_data speciesSums sampleSums nspecies species.names sampleNames sample.names getSamples getSpecies rank.names getTaxa sample.variables getVariable merge_species otuTable speciesarerows speciesAreRows plot_richness_estimates import_qiime_sampleData filterfunSample genefilterSample prune_species subset_species tipglom taxglom tre show_mothur_list_cutoffs sam_data<- sampleData<- tre<- speciesAreRows<- otuTable<- taxTab<-
#' @aliases deprecated_phyloseq_function plot_taxa_bar taxaplot taxtab taxTab sampleData samData sam_data speciesSums sampleSums nspecies species.names sampleNames sample.names getSamples getSpecies rank.names getTaxa sample.variables getVariable merge_species otuTable speciesarerows speciesAreRows plot_richness_estimates import_qiime_sampleData filterfunSample genefilterSample prune_species subset_species tipglom taxglom tre show_mothur_list_cutoffs sam_data<- sampleData<- tre<- speciesAreRows<- otuTable<- taxTab<-
#' @details
#' \tabular{rl}{
#'   \code{plot_taxa_bar} \tab now a synonym for \code{\link{plot_bar}}\cr
#'   \code{taxaplot} \tab now a synonym for \code{\link{plot_bar}}\cr
#'   \code{taxtab} \tab now a synonym for \code{\link{tax_table}}\cr
#'   \code{taxTab} \tab now a synonym for \code{\link{tax_table}}\cr
#'   \code{sampleData} \tab now a synonym for \code{\link{sample_data}}\cr
#'   \code{samData} \tab now a synonym for \code{\link{sample_data}}\cr
#'   \code{sam_data} \tab now a synonym for \code{\link{sample_data}}\cr
#'   \code{speciesSums} \tab now a synonym for \code{\link{taxa_sums}}\cr
#'   \code{sampleSums} \tab now a synonym for \code{\link{sample_sums}}\cr
#'   \code{nspecies} \tab now a synonym for \code{\link{ntaxa}}\cr
#'   \code{species.names} \tab now a synonym for \code{\link{taxa_names}}\cr
#'   \code{sampleNames} \tab now a synonym for \code{\link{sample_names}}\cr
#'   \code{sample.names} \tab now a synonym for \code{\link{sample_names}}\cr
#'   \code{getSamples} \tab now a synonym for \code{\link{get_sample}}\cr
#'   \code{getSpecies} \tab now a synonym for \code{\link{get_taxa}}\cr
#'   \code{rank.names} \tab now a synonym for \code{\link{rank_names}}\cr
#'   \code{getTaxa} \tab now a synonym for \code{\link{get_taxa_unique}}\cr
#'   \code{sample.variables} \tab now a synonym for \code{\link{sample_variables}}\cr
#'   \code{getVariable} \tab now a synonym for \code{\link{get_variable}}\cr
#'   \code{merge_species} \tab now a synonym for \code{\link{merge_taxa}}\cr
#'   \code{otuTable} \tab now a synonym for \code{\link{otu_table}}\cr
#'   \code{speciesarerows} \tab now a synonym for \code{\link{taxa_are_rows}}\cr
#'   \code{speciesAreRows} \tab now a synonym for \code{\link{taxa_are_rows}}\cr
#'   \code{plot_richness_estimates} \tab now a synonym for \code{\link{plot_richness}}\cr
#'   \code{import_qiime_sampleData} \tab now a synonym for \code{\link{import_qiime_sample_data}}\cr
#'   \code{filterfunSample} \tab now a synonym for \code{\link{filterfun_sample}}\cr
#'   \code{genefilterSample} \tab now a synonym for \code{\link{genefilter_sample}}\cr
#'   \code{prune_species} \tab now a synonym for \code{\link{prune_taxa}}\cr
#'   \code{subset_species} \tab now a synonym for \code{\link{subset_taxa}}\cr
#'   \code{tipglom} \tab now a synonym for \code{\link{tip_glom}}\cr
#'   \code{taxglom} \tab now a synonym for \code{\link{tax_glom}}\cr
#'   \code{tre} \tab now a synonym for \code{\link{phy_tree}}\cr
#'   \code{show_mothur_list_cutoffs} \tab now a synonym for \code{\link{show_mothur_cutoffs}}\cr
#'   \code{sam_data<-} \tab now a synonym for \code{\link{sample_data<-}}\cr
#'   \code{sampleData<-} \tab now a synonym for \code{\link{sample_data<-}}\cr
#'   \code{tre<-} \tab now a synonym for \code{\link{phy_tree<-}}\cr
#'   \code{speciesAreRows<-} \tab now a synonym for \code{\link{taxa_are_rows<-}}\cr
#'   \code{otuTable<-} \tab now a synonym for \code{\link{otu_table<-}}\cr
#'   \code{taxTab<-} \tab now a synonym for \code{\link{tax_table<-}}\cr
#' }
#'
deprecated_phyloseq_function <- function(x, value, ...){return(NULL)}
plot_taxa_bar <- function(...){.Deprecated("plot_bar", package="phyloseq");return(plot_bar(...))}
taxaplot <- function(...){.Deprecated("plot_bar", package="phyloseq");return(plot_bar(...))}
taxtab <- function(...){.Deprecated("tax_table", package="phyloseq");return(tax_table(...))}
taxTab <- function(...){.Deprecated("tax_table", package="phyloseq");return(tax_table(...))}
sampleData <- function(...){.Deprecated("sample_data", package="phyloseq");return(sample_data(...))}
samData <- function(...){.Deprecated("sample_data", package="phyloseq");return(sample_data(...))}
sam_data <- function(...){.Deprecated("sample_data", package="phyloseq");return(sample_data(...))}
speciesSums <- function(...){.Deprecated("taxa_sums", package="phyloseq");return(taxa_sums(...))}
sampleSums <- function(...){.Deprecated("sample_sums", package="phyloseq");return(sample_sums(...))}
nspecies <- function(...){.Deprecated("ntaxa", package="phyloseq");return(ntaxa(...))}
species.names <- function(...){.Deprecated("taxa_names", package="phyloseq");return(taxa_names(...))}
sampleNames <- function(...){.Deprecated("sample_names", package="phyloseq");return(sample_names(...))}
sample.names <- function(...){.Deprecated("sample_names", package="phyloseq");return(sample_names(...))}
getSamples <- function(...){.Deprecated("get_sample", package="phyloseq");return(get_sample(...))}
getSpecies <- function(...){.Deprecated("get_taxa", package="phyloseq");return(get_taxa(...))}
rank.names <- function(...){.Deprecated("rank_names", package="phyloseq");return(rank_names(...))}
getTaxa <- function(...){.Deprecated("get_taxa_unique", package="phyloseq");return(get_taxa_unique(...))}
sample.variables <- function(...){.Deprecated("sample_variables", package="phyloseq");return(sample_variables(...))}
getVariable <- function(...){.Deprecated("get_variable", package="phyloseq");return(get_variable(...))}
merge_species <- function(...){.Deprecated("merge_taxa", package="phyloseq");return(merge_taxa(...))}
otuTable <- function(...){.Deprecated("otu_table", package="phyloseq");return(otu_table(...))}
speciesarerows <- function(...){.Deprecated("taxa_are_rows", package="phyloseq");return(taxa_are_rows(...))}
speciesAreRows <- function(...){.Deprecated("taxa_are_rows", package="phyloseq");return(taxa_are_rows(...))}
plot_richness_estimates <- function(...){.Deprecated("plot_richness", package="phyloseq");return(plot_richness(...))}
import_qiime_sampleData <- function(...){.Deprecated("import_qiime_sample_data", package="phyloseq");return(import_qiime_sample_data(...))}
filterfunSample <- function(...){.Deprecated("filterfun_sample", package="phyloseq");return(filterfun_sample(...))}
genefilterSample <- function(...){.Deprecated("genefilter_sample", package="phyloseq");return(genefilter_sample(...))}
prune_species <- function(...){.Deprecated("prune_taxa", package="phyloseq");return(prune_taxa(...))}
subset_species <- function(...){.Deprecated("subset_taxa", package="phyloseq");return(subset_taxa(...))}
tipglom <- function(...){.Deprecated("tip_glom", package="phyloseq");return(tip_glom(...))}
taxglom <- function(...){.Deprecated("tax_glom", package="phyloseq");return(tax_glom(...))}
tre <- function(...){.Deprecated("phy_tree", package="phyloseq");return(phy_tree(...))}
show_mothur_list_cutoffs <- function(...){.Deprecated("show_mothur_cutoffs", package="phyloseq");return(show_mothur_cutoffs(...))}
originalUniFrac <- function(...){.Deprecated("fastUniFrac", package="phyloseq");return(fastUniFrac(...))}  
"sam_data<-" <- function(x, value){
  .Deprecated("sample_data<-", package="phyloseq")
  sample_data(x) <- value
  return(x)
}
"sampleData<-" <- function(x, value){
  .Deprecated("sample_data<-", package="phyloseq")
  sample_data(x) <- value
  return(x)
}
"tre<-" <- function(x, value){
  .Deprecated("phy_tree<-", package="phyloseq")
  phy_tree(x) <- value
  return(x)
}
"speciesAreRows<-" <- function(x, value){
  .Deprecated("taxa_are_rows<-", package="phyloseq")
  taxa_are_rows(x) <- value
  return(x)
}
"otuTable<-" <- function(x, value){
  .Deprecated("otu_table<-", package="phyloseq")
  otu_table(x) <- value
  return(x)
}
"taxTab<-" <- function(x, value){
  .Deprecated("tax_table<-", package="phyloseq")
  tax_table(x) <- value
  return(x)
}
################################################################################
