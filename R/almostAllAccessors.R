################################################################################
### Accessor / subset methods.
################################################################################
################################################################################
#' Retrieve reference sequences (\code{\link[Biostrings]{XStringSet}}-class) from object.
#'
#' This is the suggested method
#' for accessing
#' the phylogenetic tree, (\code{\link[Biostrings]{XStringSet}}-class)
#' from a phyloseq data object (\code{\link{phyloseq-class}}).
#' Like other accessors (see See Also, below), the default behavior of this method
#' is to stop with an
#' error if \code{physeq} is a \code{phyloseq-class} but does not 
#' contain reference sequences (the component data type you are trying to access in this case).  
#'
#' @usage refseq(physeq, errorIfNULL=TRUE)
#' 
#' @param physeq (Required). An instance of phyloseq-class
#'  that contains a phylogenetic tree. If physeq is a phylogenetic
#'  tree (a component data class), then it is returned as-is.
#'
#' @param errorIfNULL (Optional). Logical. Should the accessor stop with 
#'  an error if the slot is empty (\code{NULL})? Default \code{TRUE}.
#'
#' @return The \code{\link[ape]{phylo}}-class object contained within \code{physeq};
#'  or NULL if \code{physeq} does not have a tree.
#'  This method stops with an error in the latter NULL case be default, which
#'  can be over-ridden by changing the value of \code{errorIfNULL} to \code{FALSE}.
#'
#' @seealso \code{\link{otu_table}}, \code{\link{sample_data}}, \code{\link{tax_table}}
#'  \code{\link{phy_tree}},
#'  \code{\link{phyloseq}}, \code{\link{merge_phyloseq}}
#' 
#' @import Biostrings
#' @export
#' @rdname refseq-methods
#' @docType methods
#'
#' @examples
#'  data(GlobalPatterns)
#'  refseq(GlobalPatterns, FALSE)
setGeneric("refseq", function(physeq, errorIfNULL=TRUE) standardGeneric("refseq"))
#' @rdname refseq-methods
#' @aliases refseq,ANY-method
setMethod("refseq", "ANY", function(physeq, errorIfNULL=TRUE){
	access(physeq, "refseq", errorIfNULL)
})
# Return as-is if already a "XStringSet" object
#' @rdname refseq-methods
#' @aliases refseq,XStringSet-method
setMethod("refseq", "XStringSet", function(physeq){ return(physeq) })
################################################################################
#' Retrieve phylogenetic tree (\code{\link[ape]{phylo}}-class) from object.
#'
#' This is the suggested method
#' for accessing
#' the phylogenetic tree, (\code{\link[ape]{phylo}}-class) from a \code{\link{phyloseq-class}}.
#' Like other accessors (see See Also, below), the default behavior of this method
#' is to stop with an
#' error if \code{physeq} is a \code{phyloseq-class} but does not 
#' contain a phylogenetic tree (the component data you are trying to access in this case).  
#' 
#' Note that the tip labels should be named to match the
#' \code{taxa_names} of the other objects to which it is going to be paired.
#' The \code{\link{phyloseq}} constructor automatically checks for
#' exact agreement in the set of species described by the phlyogenetic tree 
#' and the other components (taxonomyTable, otu_table),
#' and trims as-needed. Thus, the tip.labels in a phylo object
#' must be named to match the results of
#' \code{\link{taxa_names}} of the other objects to which it will ultimately be paired.
#'
#' @usage phy_tree(physeq, errorIfNULL=TRUE)
#' 
#' @param physeq (Required). An instance of phyloseq-class
#'  that contains a phylogenetic tree. If physeq is a phylogenetic
#'  tree (a component data class), then it is returned as-is.
#'
#' @param errorIfNULL (Optional). Logical. Should the accessor stop with 
#'  an error if the slot is empty (\code{NULL})? Default \code{TRUE}.
#'
#' @return The \code{\link[ape]{phylo}}-class object contained within \code{physeq};
#'  or NULL if \code{physeq} does not have a tree.
#'  This method stops with an error in the latter NULL case be default, which
#'  can be over-ridden by changing the value of \code{errorIfNULL} to \code{FALSE}.
#'
#' @seealso \code{\link{otu_table}}, \code{\link{sample_data}}, \code{\link{tax_table}}
#'  \code{\link{refseq}},
#'  \code{\link{phyloseq}}, \code{\link{merge_phyloseq}}
#' 
#' @export
#' @rdname phy_tree-methods
#' @docType methods
#'
#' @examples
#'  data(GlobalPatterns)
#'  phy_tree(GlobalPatterns)
setGeneric("phy_tree", function(physeq, errorIfNULL=TRUE) standardGeneric("phy_tree"))
#' @rdname phy_tree-methods
#' @aliases phy_tree,ANY-method
setMethod("phy_tree", "ANY", function(physeq, errorIfNULL=TRUE){
	access(physeq, "phy_tree", errorIfNULL)
})
# Return as-is if already a "phylo" object
#' @rdname phy_tree-methods
#' @aliases phy_tree,phylo-method
setMethod("phy_tree", "phylo", function(physeq){ return(physeq) })
################################################################################
#' Access taxa_are_rows slot from otu_table objects.
#'
#' @usage taxa_are_rows(physeq)
#'
#' @param physeq (Required). \code{\link{phyloseq-class}}, or \code{\link{otu_table-class}}.
#'
#' @return A logical indicating the orientation of the otu_table.
#'
#' @seealso \code{\link{otu_table}}
#' @rdname taxa_are_rows-methods
#' @docType methods
#' @export
#' @aliases taxa_are_rows taxa_are_rows
setGeneric("taxa_are_rows", function(physeq) standardGeneric("taxa_are_rows"))
#' @rdname taxa_are_rows-methods
#' @aliases taxa_are_rows,ANY-method
setMethod("taxa_are_rows", "ANY", function(physeq){NULL})
#' @rdname taxa_are_rows-methods
#' @aliases taxa_are_rows,otu_table-method
setMethod("taxa_are_rows", "otu_table", function(physeq){physeq@taxa_are_rows})
#' @rdname taxa_are_rows-methods
#' @aliases taxa_are_rows,phyloseq-method
setMethod("taxa_are_rows", "phyloseq", function(physeq){
	taxa_are_rows(otu_table(physeq))
})
################################################################################
#' Get the number of taxa/species.
#'
#' @usage ntaxa(physeq)
#'
#' @param physeq \code{\link{phyloseq-class}}, \code{\link{otu_table-class}},
#'  \code{\link{taxonomyTable-class}}, or
#'  \code{\link[ape]{phylo}}
#'
#' @return An integer indicating the number of taxa / species.
#'
#' @seealso taxa_names
#'
#' @rdname ntaxa-methods
#' @docType methods
#' @export
#'
#' @examples #
#' # # From "picante" package
#' # data("phylocom")
#' # tree <- phylocom$phylo
#' # ntaxa(tree)
setGeneric("ntaxa", function(physeq) standardGeneric("ntaxa"))
#' @rdname ntaxa-methods
#' @aliases ntaxa,ANY-method
setMethod("ntaxa", "ANY", function(physeq){ return(NULL) })
#' @rdname ntaxa-methods
#' @aliases ntaxa,phyloseq-method
setMethod("ntaxa", "phyloseq", function(physeq){
	length(taxa_names(physeq))
})
#' @rdname ntaxa-methods
#' @aliases ntaxa,otu_table-method
setMethod("ntaxa", "otu_table", function(physeq){
	if( taxa_are_rows(physeq) ){
		return( length(rownames(physeq)) )
	} else {
		return( length(colnames(physeq)) )
	}
})
#' @rdname ntaxa-methods
#' @aliases ntaxa,taxonomyTable-method
setMethod("ntaxa", "taxonomyTable", function(physeq){ nrow(physeq) })
#' @rdname ntaxa-methods
#' @aliases ntaxa,phylo-method
setMethod("ntaxa", "phylo", function(physeq) length(physeq$tip.label) )
#' @rdname ntaxa-methods
#' @aliases ntaxa,XStringSet-method
setMethod("ntaxa", "XStringSet", function(physeq) length(physeq) )
################################################################################
#' Get species / taxa names.
#'
#' @usage taxa_names(physeq)
#'
#' @param physeq \code{\link{phyloseq-class}}, \code{\link{otu_table-class}},
#'  \code{\link{taxonomyTable-class}}, or
#'  \code{\link[ape]{phylo}}
#'
#' @return A character vector of the names of the species in \code{physeq}.
#'
#' @seealso ntaxa
#'
#' @rdname taxa_names-methods
#' @docType methods
#' @export
#'
#' @examples #
#' # # From "picante" package
#' # data("phylocom")
#' # tree <- phylocom$phylo
#' # OTU1 <- otu_table(phylocom$sample, taxa_are_rows=FALSE)
#' # taxa_names(tree)
#' # taxa_names(OTU1)
#' # physeq1 <- phyloseq(OTU1, tree)
#' # taxa_names(physeq1)
setGeneric("taxa_names", function(physeq) standardGeneric("taxa_names"))	
#' @rdname taxa_names-methods
#' @aliases taxa_names,ANY-method
setMethod("taxa_names", "ANY", function(physeq){ return(NULL) })
#' @rdname taxa_names-methods
#' @aliases taxa_names,phyloseq-method
setMethod("taxa_names", "phyloseq", function(physeq){
	# Return the union of all species in the components
	unique(unlist(lapply(splat.phyloseq.objects(physeq), taxa_names)))	
})
#' @rdname taxa_names-methods
#' @aliases taxa_names,otu_table-method
setMethod("taxa_names", "otu_table", function(physeq){
	if( taxa_are_rows(physeq) ){
		return( rownames(physeq) )
	} else {
		return( colnames(physeq) )
	}
})
#' @rdname taxa_names-methods
#' @aliases taxa_names,taxonomyTable-method
setMethod("taxa_names", "taxonomyTable", function(physeq) rownames(physeq) )
#' @rdname taxa_names-methods
#' @aliases taxa_names,sample_data-method
setMethod("taxa_names", "sample_data", function(physeq) NULL )
#' @rdname taxa_names-methods
#' @aliases taxa_names,phylo-method
setMethod("taxa_names", "phylo", function(physeq) physeq$tip.label )
#' @rdname taxa_names-methods
#' @aliases taxa_names,XStringSet-method
setMethod("taxa_names", "XStringSet", function(physeq) names(physeq) )
################################################################################
#' Get the number of samples.
#'
#' @usage nsamples(physeq)
#'
#' @param physeq A \code{\link{phyloseq-class}}, \code{\link{sample_data}},
#'  or \code{\link{otu_table-class}}.
#'
#' @return An integer indicating the total number of samples.
#'
#' @seealso \code{\link{taxa_names}}, \code{\link{sample_names}},
#'  \code{\link{ntaxa}}
#'
#' @rdname nsamples-methods
#' @docType methods
#' @export
#'
#' @examples #
#' # # From "picante" package
#' # data("phylocom")
#' # tree <- phylocom$phylo
#' # OTU1 <- otu_table(phylocom$sample, taxa_are_rows=FALSE)
#' # nsamples(OTU1)
#' # physeq1 <- phyloseq(OTU1, tree)
#' # nsamples(physeq1)
setGeneric("nsamples", function(physeq) standardGeneric("nsamples"))
#' @rdname nsamples-methods
#' @aliases nsamples,ANY-method
setMethod("nsamples", "ANY", function(physeq){ return(NULL) })
#' @rdname nsamples-methods
#' @aliases nsamples,phyloseq-method
setMethod("nsamples", "phyloseq", function(physeq){
	length(sample_names(physeq))
})
#' @rdname nsamples-methods
#' @aliases nsamples,otu_table-method
setMethod("nsamples", "otu_table", function(physeq){
	if( taxa_are_rows(physeq) ){
		return( length(colnames(physeq)) )
	} else {
		return( length(rownames(physeq)) )
	}	
})
#' @rdname nsamples-methods
#' @aliases nsamples,sample_data-method
setMethod("nsamples", "sample_data", function(physeq) nrow(physeq) )
################################################################################
#' Get sample names.
#'
#' @usage sample_names(physeq)
#'
#' @param physeq (Required). A \code{\link{phyloseq-class}}, \code{\link{sample_data}},
#'  or \code{\link{otu_table-class}}.
#'
#' @return A character vector. The names of the samples in \code{physeq}.
#'
#' @seealso \code{\link{taxa_names}}, \code{\link{nsamples}}
#' 
#' @aliases sample_names sampleNames
#'
#' @rdname sample_names-methods
#' @docType methods
#' @export
#'
#' @examples #
#' # # From "picante" package
#' # data(GlobalPatterns)
#' # sample_names(GlobalPatterns)
setGeneric("sample_names", function(physeq) standardGeneric("sample_names"))
# Unless otherwise specified, this should return a value of NULL
# That way, objects that do not explicitly describe samples all
# behave in the same (returning NULL) way.
#' @rdname sample_names-methods
#' @aliases sample_names,ANY-method
setMethod("sample_names", "ANY", function(physeq){ return(NULL) })
#' @rdname sample_names-methods
#' @aliases sample_names,phyloseq-method
setMethod("sample_names", "phyloseq", function(physeq){
	# Return the union of all samples in the components
	unique(unlist(lapply(splat.phyloseq.objects(physeq), sample_names)))
})
#' @rdname sample_names-methods
#' @aliases sample_names,sample_data-method
setMethod("sample_names", "sample_data", function(physeq) rownames(physeq) )
#' @rdname sample_names-methods
#' @aliases sample_names,otu_table-method
setMethod("sample_names", "otu_table", function(physeq){
	if( taxa_are_rows(physeq) ){
		return( colnames(physeq) )
	} else {
		return( rownames(physeq) )
	}
})
################################################################################
#' Returns all abundance values for species \code{i}.
#'
#' This is a simple accessor function for investigating 
#' a single species-of-interest. 
#'
#' @usage get_sample(physeq, i)
#' @param physeq (Required). \code{\link{otu_table-class}}, or \code{\link{phyloseq-class}}.
#' @param i (Required). A single taxa/species/OTU ID for which you want
#'  to know the abundance in each sample.
#'
#' @return An integer vector of the abundance values for 
#' each sample in \code{physeq} for species \code{i}
#' 
#' @seealso
#'  \code{\link{get_taxa}}
#'  \code{\link{taxa_names}}
#'  \code{\link{sample_names}}
#'
#' @rdname get_sample-methods
#' @docType methods
#' @export
#'
#' @examples
#' data(esophagus)
#' taxa_names(esophagus)
#' get_sample(esophagus, "59_5_19")
setGeneric("get_sample", function(physeq, i) standardGeneric("get_sample"))
################################################################################
#' @aliases get_sample,otu_table-method
#' @rdname get_sample-methods
setMethod("get_sample", "otu_table", function(physeq, i){
	if( taxa_are_rows(physeq) ){
		as(physeq, "matrix")[i, ]
	} else {
		as(physeq, "matrix")[, i]
	}
})
################################################################################
#' @aliases get_sample,phyloseq-method
#' @rdname get_sample-methods
setMethod("get_sample", "phyloseq", function(physeq, i){
	get_sample(otu_table(physeq), i)
})
################################################################################
#' Returns all abundance values of sample \code{i}.
#'
#' This is a simple accessor function for investigating 
#' a single sample-of-interest. 
#'
#' @usage get_taxa(physeq, i)
#' @param physeq (Required). \code{\link{otu_table-class}}, or \code{\link{phyloseq-class}}.
#' @param i (Required). A single sample for which you want
#'  to know the abundance of each species. Can be integer
#'  for index value, or sample name.
#'
#' @return An integer vector of the abundance values for 
#' each species in \code{physeq} for sample \code{i}
#' 
#' @seealso
#'  \code{\link{get_sample}}
#'  \code{\link{taxa_names}}
#'  \code{\link{sample_names}}
#'
#' @rdname get_taxa-methods
#' @docType methods
#' @export
#'
#' @examples
#' data(esophagus)
#' sample_names(esophagus)
#' get_taxa(esophagus, "B")
setGeneric("get_taxa", function(physeq, i) standardGeneric("get_taxa"))
#' @aliases get_taxa,otu_table-method
#' @rdname get_taxa-methods
setMethod("get_taxa", "otu_table", function(physeq, i){
	if( taxa_are_rows(physeq) ){
		as(physeq, "matrix")[, i]
	} else {
		as(physeq, "matrix")[i, ]
	}
})
#' @aliases get_taxa,phyloseq-method
#' @rdname get_taxa-methods
setMethod("get_taxa", "phyloseq", function(physeq, i){
	get_taxa(otu_table(physeq), i)
})
################################################################################
#' Retrieve the names of the taxonomic ranks 
#'
#' This is a simple accessor function to make it more convenient to determine
#' the taxonomic ranks that are available in a given \code{\link{phyloseq-class}}
#' object. 
#'
#' @usage rank_names(physeq, errorIfNULL=TRUE)
#' 
#' @param physeq (Required). 
#'  \code{\link{taxonomyTable-class}}, or \code{\link{phyloseq-class}}.
#'
#' @param errorIfNULL (Optional). Logical. Should the accessor stop with 
#'  an error if the slot is empty (\code{NULL})? Default \code{TRUE}.
#'
#' @return Character vector. The names of the available taxonomic ranks.
#' 
#' @seealso
#'  \code{\link{get_taxa}}
#'  \code{\link{taxa_names}}
#'  \code{\link{sample_names}}
#'
#' @export
#'
#' @examples
#' data(enterotype)
#' rank_names(enterotype)
rank_names <- function(physeq, errorIfNULL=TRUE){
	colnames(tax_table(physeq, errorIfNULL))	
}
################################################################################
#' Get a unique vector of the observed taxa at a particular taxonomic rank
#'
#' This is a simple accessor function to make it more convenient to determine
#' the different taxa present for a particular taxonomic rank
#' in a given \code{\link{phyloseq-class}} object. 
#'
#' @usage get_taxa_unique(physeq, taxonomic.rank=rank_names(physeq)[1], errorIfNULL=TRUE)
#' 
#' @param physeq (Required). \code{\link{taxonomyTable-class}}, or \code{\link{phyloseq-class}}.
#'
#' @param taxonomic.rank (Optional). Character. The taxonomic rank to use. Must select
#'  from the set indicated by \code{get_taxa_unique}. Default is
#'  to take the first column of the \code{taxonomyTable} component.
#'
#' @param errorIfNULL (Optional). Logical. Should the accessor stop with 
#'  an error if the slot is empty (\code{NULL})? Default \code{TRUE}.
#'
#' @return Character vector. Unique vector of the observed taxa 
#'  at a particular taxonomic rank
#' 
#' @seealso
#'  \code{\link{get_taxa}}
#'  \code{\link{taxa_names}}
#'  \code{\link{sample_names}}
#'
#' @export
#'
#' @examples
#' data(enterotype)
#' get_taxa_unique(enterotype)
#' data(GlobalPatterns)
#' get_taxa_unique(GlobalPatterns, "Family")
get_taxa_unique <- function(physeq, taxonomic.rank=rank_names(physeq)[1], errorIfNULL=TRUE){
	unique(as(tax_table(physeq, errorIfNULL)[, taxonomic.rank], "character"))
}
################################################################################
#' Get the sample variables present in sample_data
#'
#' This is a simple accessor function to make it more convenient to determine
#' the sample variable names of a particular \code{\link{phyloseq-class}} object. 
#'
#' @usage sample_variables(physeq, errorIfNULL=TRUE)
#' 
#' @param physeq (Required). \code{\link{sample_data-class}}, or \code{\link{phyloseq-class}}.
#'
#' @param errorIfNULL (Optional). Logical. Should the accessor stop with 
#'  an error if the slot is empty (\code{NULL})? Default \code{TRUE}.
#'
#' @return Character vector. The names of the variables in the sample_data
#'  data.frame. Essentially the column names. Useful for selecting model 
#'  and graphics parameters that interact with sample_data.
#' 
#' @seealso
#'  \code{\link{get_taxa}}
#'  \code{\link{taxa_names}}
#'  \code{\link{sample_names}}
#'
#' @export
#'
#' @examples
#' data(enterotype)
#' sample_variables(enterotype)
sample_variables <- function(physeq, errorIfNULL=TRUE){
	colnames(sample_data(physeq, errorIfNULL))
}
################################################################################
#' Get the values for a particular variable in sample_data
#'
#' This is a simple accessor function for streamlining access
#' to values/vectors/factors/etc contained in the sample_data.
#'
#' @usage get_variable(physeq, varName)
#' 
#' @param physeq (Required). \code{\link{sample_data-class}}, or \code{\link{phyloseq-class}}.
#'
#' @param varName (Required). Character string of the variable name in \code{sample_data}.
#'  Use \code{sample_variables(physeq)} for available variables in your object.
#'
#' @return Data. The clas of the data depends on what the contents of sample_data.
#' 
#' @seealso
#'  \code{\link{get_taxa}}
#'  \code{\link{taxa_names}}
#'  \code{\link{sample_names}}
#'
#'  \code{\link{sample_variables}}
#'
#' @export
#'
#' @examples
#' # Load the GlobalPatterns dataset into the workspace environment
#' data(GlobalPatterns)
#' # Look at the different values for SampleType 
#' get_variable(GlobalPatterns, "SampleType")
get_variable <- function(physeq, varName){
	if( is.null(sample_data(physeq, FALSE)) ){
		stop("Your phyloseq data object does not have a sample-data component\n",
			"Try ?sample_data for more details.")
	}
	return( as(sample_data(physeq), "data.frame")[, varName] )
}
################################################################################
