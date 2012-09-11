# Deprecated function definitions and roxygen headers (source).
# Backward compatibility aliases



#' @rdname plot_richness
#' @aliases plot_richness
#' @export
plot_richness_estimates <- plot_richness


#' @rdname import_qiime_sample_data
#' @aliases import_qiime_sample_data
#' @export
import_qiime_sampleData <- import_qiime_sample_data

#' @rdname filterfun_sample
#' @aliases filterfun_sample
#' @export
filterfunSample <- filterfun_sample



#' @rdname genefilter_sample-methods
#' @aliases genefilter_sample
#' @export
genefilterSample <- genefilter_sample



#' @rdname prune_taxa-methods
#' @aliases prune_taxa
#' @export
prune_species <- prune_taxa



#' @rdname subset_taxa-methods
#' @docType methods
#' @export
subset_species <- subset_taxa



#' @rdname tip_glom-methods
#' @aliases tip_glom
#' @export
tipglom <- tip_glom



#' @rdname tax_glom
#' @aliases tax_glom
#' @export
taxglom <- tax_glom



#' @rdname phy_tree-methods
#' @aliases phy_tree
#' @export
tre <- phy_tree

#' @rdname assign-phy_tree
#' @aliases assign-phy_tree phy_tree<-
#' @export
"tre<-" <- function(x, value){
	phy_tree(x) <- value
	return(x)
}


#' @rdname taxa_are_rows-methods
#' @aliases taxa_are_rows
#' @export
speciesarerows <- taxa_are_rows

#' @rdname taxa_are_rows-methods
#' @aliases taxa_are_rows
#' @export
speciesAreRows <- taxa_are_rows

#' @rdname assign-taxa_are_rows
#' @aliases assign-taxa_are_rows
#' @export
"speciesAreRows<-" <- function(x, value){
	taxa_are_rows(x) <- value
	return(x)
}


#' @rdname ntaxa-methods
#' @aliases ntaxa
#' @export
nspecies <- ntaxa



#' @rdname taxa_names-methods
#' @aliases taxa_names
#' @export
species.names <- taxa_names



#' @rdname sample_names-methods
#' @aliases sample_names
#' @export
sampleNames <- sample_names

#' @rdname sample_names-methods
#' @aliases sample_names
#' @export
sample.names <- sample_names



#' @rdname get_sample-methods
#' @aliases get_sample
#' @export
getSamples <- get_sample



#' @rdname get_taxa-methods
#' @aliases get_taxa
#' @export
getSpecies <- get_taxa



#' @rdname rank_names
#' @aliases rank_names
#' @export
rank.names <- rank_names



#' @rdname get_taxa_unique
#' @aliases get_taxa_unique
#' @export
getTaxa <- get_taxa_unique



#' @rdname sample_variables
#' @aliases sample_variables
#' @export
sample.variables <- sample_variables



#' @rdname get_variable
#' @aliases get_variable
#' @export
getVariable <- get_variable



#' @rdname merge_taxa-methods
#' @aliases merge_taxa
#' @export
merge_species <- merge_taxa



#' @rdname otu_table-methods
#' @aliases otu_table
#' @export
otuTable <- otu_table

#' @rdname assign-otu_table
#' @aliases assign-otu_table otu_table<- otuTable<-
#' @export
"otuTable<-" <- function(x, value){
	otu_table(x) <- value
	return(x)
}



#' @rdname taxa_sums
#' @aliases taxa_sums
#' @export
speciesSums <- taxa_sums




#' @rdname sample_sums
#' @aliases sample_sums
#' @export
sampleSums <- sample_sums



#' @rdname sample_data-methods
#' @aliases sample_data
#' @export
sampleData <- sample_data

#' @rdname sample_data-methods
#' @aliases sample_data
#' @export
samData <- sample_data

#' @rdname assign-sample_data
#' @aliases sample_data<- 
#' @export
"sampleData<-" <- function(x, value){
	sample_data(x) <- value
	return(x)
}



#' @rdname tax_table-methods
#' @aliases tax_table
#' @export
taxtab <- tax_table

#' @rdname tax_table-methods
#' @aliases tax_table
#' @export
taxTab <- tax_table

#' @export
#' @rdname assign-tax_table
#' @aliases tax_table<-
"taxTab<-" <- function(x, value){
	tax_table(x) <- value
	return(x)
}


