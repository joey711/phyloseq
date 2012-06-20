################################################################################
### Accessor / subset methods.
################################################################################
#' Get phylogenetic tree from object.
#'
#' This is the main method suggested 
#' for accessing
#' the phylogenetic tree, (\code{\link[ape]{phylo}}-class) from a \code{\link{phyloseq-class}}.
#' Like other accessors (see See Also, below), the default behavior of this method
#' is to stop with an
#' error if \code{physeq} is a \code{phyloseq-class} but does not 
#' contain a phylogenetic tree.  
#' 
#' Note that the tip labels should be named to match the
#' \code{species.names} of the other objects to which it is going to be paired.
#' The \code{\link{phyloseq}} constructor automatically checks for
#' exact agreement in the set of species described by the phlyogenetic tree 
#' and the other components (taxonomyTable, otuTable),
#' and trims as-needed. Thus, the tip.labels in a phylo object
#' must be named to match the results of
#' \code{\link{species.names}} of the other objects to which it will ultimately be paired.
#'
#' @usage tre(physeq, errorIfNULL=TRUE)
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
#' @seealso \code{\link{otuTable}}, \code{\link{sampleData}}, \code{\link{taxTab}}
#'  \code{\link{phyloseq}}, \code{\link{merge_phyloseq}}
#' 
#' @export
#' @rdname tre-methods
#' @docType methods
#'
#' @examples
#' # data(GlobalPatterns)
#' # tre(GlobalPatterns)
setGeneric("tre", function(physeq, errorIfNULL=TRUE) standardGeneric("tre"))
#' @rdname tre-methods
#' @aliases tre,ANY-method
setMethod("tre", "ANY", function(physeq, errorIfNULL=TRUE){
	access(physeq, "tre", errorIfNULL)
})
# Return as-is if already a "phylo" object
#' @rdname tre-methods
#' @aliases tre,phylo-method
setMethod("tre", "phylo", function(physeq){ return(physeq) })
################################################################################
#' Access speciesAreRows slot from otuTable objects.
#'
#' @usage speciesarerows(physeq)
#'
#' @param physeq (Required). \code{\link{phyloseq-class}}, or \code{\link{otuTable-class}}.
#'
#' @return A logical indicating the orientation of the otuTable.
#'
#' @seealso \code{\link{otuTable}}
#' @rdname speciesAreRows-methods
#' @docType methods
#' @export
#' @aliases speciesAreRows speciesarerows
setGeneric("speciesAreRows", function(physeq) standardGeneric("speciesAreRows"))
#' @rdname speciesAreRows-methods
#' @aliases speciesAreRows,ANY-method
setMethod("speciesAreRows", "ANY", function(physeq){NULL})
#' @rdname speciesAreRows-methods
#' @aliases speciesAreRows,otuTable-method
setMethod("speciesAreRows", "otuTable", function(physeq){physeq@speciesAreRows})
#' @rdname speciesAreRows-methods
#' @aliases speciesAreRows,phyloseq-method
setMethod("speciesAreRows", "phyloseq", function(physeq){
	speciesAreRows(otuTable(physeq))
})
#' @aliases speciesAreRows speciesarerows
#' @export
speciesarerows <- speciesAreRows
################################################################################
#' Get the number of taxa/species.
#'
#' @usage nspecies(physeq)
#'
#' @param physeq \code{\link{phyloseq-class}}, \code{\link{otuTable-class}},
#'  \code{\link{taxonomyTable-class}}, or
#'  \code{\link[ape]{phylo}}
#'
#' @return An integer indicating the number of taxa / species.
#'
#' @seealso species.names
#'
#' @rdname nspecies-methods
#' @docType methods
#' @export
#'
#' @examples #
#' # # From "picante" package
#' # data("phylocom")
#' # tree <- phylocom$phylo
#' # nspecies(tree)
setGeneric("nspecies", function(physeq) standardGeneric("nspecies"))
#' @rdname nspecies-methods
#' @aliases nspecies,ANY-method
setMethod("nspecies", "ANY", function(physeq){ return(NULL) })
#' @rdname nspecies-methods
#' @aliases nspecies,phyloseq-method
setMethod("nspecies", "phyloseq", function(physeq){
	length(species.names(physeq))
})
#' @rdname nspecies-methods
#' @aliases nspecies,otuTable-method
setMethod("nspecies", "otuTable", function(physeq){
	if( speciesAreRows(physeq) ){
		return( length(rownames(physeq)) )
	} else {
		return( length(colnames(physeq)) )
	}
})
#' @rdname nspecies-methods
#' @aliases nspecies,taxonomyTable-method
setMethod("nspecies", "taxonomyTable", function(physeq){ nrow(physeq) })
#' @rdname nspecies-methods
#' @aliases nspecies,phylo-method
setMethod("nspecies", "phylo", function(physeq) length(physeq$tip.label) )
################################################################################
#' Get species / taxa names.
#'
#' @usage species.names(physeq)
#'
#' @param physeq \code{\link{phyloseq-class}}, \code{\link{otuTable-class}},
#'  \code{\link{taxonomyTable-class}}, or
#'  \code{\link[ape]{phylo}}
#'
#' @return A character vector of the names of the species in \code{physeq}.
#'
#' @seealso nspecies
#'
#' @rdname species.names-methods
#' @docType methods
#' @export
#'
#' @examples #
#' # # From "picante" package
#' # data("phylocom")
#' # tree <- phylocom$phylo
#' # OTU1 <- otuTable(phylocom$sample, speciesAreRows=FALSE)
#' # species.names(tree)
#' # species.names(OTU1)
#' # physeq1 <- phyloseq(OTU1, tree)
#' # species.names(physeq1)
setGeneric("species.names", function(physeq) standardGeneric("species.names"))	
#' @rdname species.names-methods
#' @aliases species.names,ANY-method
setMethod("species.names", "ANY", function(physeq){ return(NULL) })
#' @rdname species.names-methods
#' @aliases species.names,phyloseq-method
setMethod("species.names", "phyloseq", function(physeq){
	# Return the union of all species in the components
	unique(unlist(lapply(splat.phyloseq.objects(physeq), species.names)))	
})
#' @rdname species.names-methods
#' @aliases species.names,otuTable-method
setMethod("species.names", "otuTable", function(physeq){
	if( speciesAreRows(physeq) ){
		return( rownames(physeq) )
	} else {
		return( colnames(physeq) )
	}
})
#' @rdname species.names-methods
#' @aliases species.names,taxonomyTable-method
setMethod("species.names", "taxonomyTable", function(physeq) rownames(physeq) )
#' @rdname species.names-methods
#' @aliases species.names,sampleData-method
setMethod("species.names", "sampleData", function(physeq) NULL )
#' @rdname species.names-methods
#' @aliases species.names,phylo-method
setMethod("species.names", "phylo", function(physeq) physeq$tip.label )
################################################################################
#' Get the number of samples.
#'
#' @usage nsamples(physeq)
#'
#' @param physeq A \code{\link{phyloseq-class}}, \code{\link{sampleData}},
#'  or \code{\link{otuTable-class}}.
#'
#' @return An integer indicating the total number of samples.
#'
#' @seealso \code{\link{species.names}}, \code{\link{sample.names}},
#'  \code{\link{nspecies}}
#'
#' @rdname nsamples-methods
#' @docType methods
#' @export
#'
#' @examples #
#' # # From "picante" package
#' # data("phylocom")
#' # tree <- phylocom$phylo
#' # OTU1 <- otuTable(phylocom$sample, speciesAreRows=FALSE)
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
	length(sample.names(physeq))
})
#' @rdname nsamples-methods
#' @aliases nsamples,otuTable-method
setMethod("nsamples", "otuTable", function(physeq){
	if( speciesAreRows(physeq) ){
		return( length(colnames(physeq)) )
	} else {
		return( length(rownames(physeq)) )
	}	
})
#' @rdname nsamples-methods
#' @aliases nsamples,sampleData-method
setMethod("nsamples", "sampleData", function(physeq) nrow(physeq) )
################################################################################
#' Get sample names.
#'
#' @usage sample.names(physeq)
#'
#' @param physeq (Required). A \code{\link{phyloseq-class}}, \code{\link{sampleData}},
#'  or \code{\link{otuTable-class}}.
#'
#' @return A character vector. The names of the samples in \code{physeq}.
#'
#' @seealso \code{\link{species.names}}, \code{\link{nsamples}}
#' 
#' @aliases sample.names sampleNames
#'
#' @rdname sample.names-methods
#' @docType methods
#' @export
#'
#' @examples #
#' # # From "picante" package
#' # data(GlobalPatterns)
#' # sample.names(GlobalPatterns)
setGeneric("sample.names", function(physeq) standardGeneric("sample.names"))
# Unless otherwise specified, this should return a value of NULL
# That way, objects that do not explicitly describe samples all
# behave in the same (returning NULL) way.
#' @rdname sample.names-methods
#' @aliases sample.names,ANY-method
setMethod("sample.names", "ANY", function(physeq){ return(NULL) })
#' @rdname sample.names-methods
#' @aliases sample.names,phyloseq-method
setMethod("sample.names", "phyloseq", function(physeq){
	# Return the union of all samples in the components
	unique(unlist(lapply(splat.phyloseq.objects(physeq), sample.names)))
})
#' @rdname sample.names-methods
#' @aliases sample.names,sampleData-method
setMethod("sample.names", "sampleData", function(physeq) rownames(physeq) )
#' @rdname sample.names-methods
#' @aliases sample.names,otuTable-method
setMethod("sample.names", "otuTable", function(physeq){
	if( speciesAreRows(physeq) ){
		return( colnames(physeq) )
	} else {
		return( rownames(physeq) )
	}
})
#' @aliases sample.names
#' @export
sampleNames <- sample.names
################################################################################
#' Returns all abundance values for species \code{i}.
#'
#' This is a simple accessor function for investigating 
#' a single species-of-interest. 
#'
#' @usage getSamples(physeq, i)
#' @param physeq (Required). \code{\link{otuTable-class}}, or \code{\link{phyloseq-class}}.
#' @param i (Required). A single taxa/species/OTU ID for which you want
#'  to know the abundance in each sample.
#'
#' @return An integer vector of the abundance values for 
#' each sample in \code{physeq} for species \code{i}
#' 
#' @seealso getSpecies species.names sample.names
#'
#' @rdname getSamples-methods
#' @docType methods
#' @export
#'
#' @examples
#' data(esophagus)
#' species.names(esophagus)
#' getSamples(esophagus, "59_5_19")
setGeneric("getSamples", function(physeq, i) standardGeneric("getSamples"))
################################################################################
#' @aliases getSamples,otuTable-method
#' @rdname getSamples-methods
setMethod("getSamples", "otuTable", function(physeq, i){
	if( speciesAreRows(physeq) ){
		as(physeq, "matrix")[i, ]
	} else {
		as(physeq, "matrix")[, i]
	}
})
################################################################################
#' @aliases getSamples,phyloseq-method
#' @rdname getSamples-methods
setMethod("getSamples", "phyloseq", function(physeq, i){
	getSamples(otuTable(physeq), i)
})
################################################################################
#' Returns all abundance values of sample \code{i}.
#'
#' This is a simple accessor function for investigating 
#' a single sample-of-interest. 
#'
#' @usage getSpecies(physeq, i)
#' @param physeq (Required). \code{\link{otuTable-class}}, or \code{\link{phyloseq-class}}.
#' @param i (Required). A single sample for which you want
#'  to know the abundance of each species. Can be integer
#'  for index value, or sample name.
#'
#' @return An integer vector of the abundance values for 
#' each species in \code{physeq} for sample \code{i}
#' 
#' @seealso getSpecies species.names sample.names
#'
#' @rdname getSpecies-methods
#' @docType methods
#' @export
#'
#' @examples
#' data(esophagus)
#' sample.names(esophagus)
#' getSpecies(esophagus, "B")
setGeneric("getSpecies", function(physeq, i) standardGeneric("getSpecies"))
#' @aliases getSpecies,otuTable-method
#' @rdname getSpecies-methods
setMethod("getSpecies", "otuTable", function(physeq, i){
	if( speciesAreRows(physeq) ){
		as(physeq, "matrix")[, i]
	} else {
		as(physeq, "matrix")[i, ]
	}
})
#' @aliases getSpecies,phyloseq-method
#' @rdname getSpecies-methods
setMethod("getSpecies", "phyloseq", function(physeq, i){
	getSpecies(otuTable(physeq), i)
})
################################################################################
#' Get the names of the taxonomic ranks 
#'
#' This is a simple accessor function to make it more convenient to determine
#' the taxonomic ranks that are available in a given \code{\link{phyloseq-class}}
#' object. 
#'
#' @usage rank.names(physeq, errorIfNULL=TRUE)
#' 
#' @param physeq (Required). 
#'  \code{\link{taxonomyTable-class}}, or \code{\link{phyloseq-class}}.
#'
#' @param errorIfNULL (Optional). Logical. Should the accessor stop with 
#'  an error if the slot is empty (\code{NULL})? Default \code{TRUE}.
#'
#' @return Character vector. The names of the available taxonomic ranks.
#' 
#' @seealso getSpecies species.names sample.names getTaxa
#'
#' @export
#'
#' @examples
#' data(enterotype)
#' rank.names(enterotype)
rank.names <- function(physeq, errorIfNULL=TRUE){
	colnames(taxTab(physeq, errorIfNULL))	
}
################################################################################
#' Get a unique vector of the observed taxa at a particular taxonomic rank
#'
#' This is a simple accessor function to make it more convenient to determine
#' the different taxa present for a particular taxonomic rank
#' in a given \code{\link{phyloseq-class}} object. 
#'
#' @usage getTaxa(physeq, taxonomic.rank=rank.names(physeq)[1], errorIfNULL=TRUE)
#' 
#' @param physeq (Required). \code{\link{taxonomyTable-class}}, or \code{\link{phyloseq-class}}.
#'
#' @param taxonomic.rank (Optional). Character. The taxonomic rank to use. Must select
#'  from the set indicated by \code{getTaxa}. Default is
#'  to take the first column of the \code{taxonomyTable} component.
#'
#' @param errorIfNULL (Optional). Logical. Should the accessor stop with 
#'  an error if the slot is empty (\code{NULL})? Default \code{TRUE}.
#'
#' @return Character vector. Unique vector of the observed taxa 
#'  at a particular taxonomic rank
#' 
#' @seealso getSpecies species.names sample.names getTaxa
#'
#' @export
#'
#' @examples
#' data(enterotype)
#' getTaxa(enterotype)
#' data(GlobalPatterns)
#' getTaxa(GlobalPatterns, "Family")
getTaxa <- function(physeq, taxonomic.rank=rank.names(physeq)[1], errorIfNULL=TRUE){
	unique(as(taxTab(physeq, errorIfNULL)[, taxonomic.rank], "character"))
}
################################################################################
#' Get the sample variables present in sampleData
#'
#' This is a simple accessor function to make it more convenient to determine
#' the sample variable names of a particular \code{\link{phyloseq-class}} object. 
#'
#' @usage sample.variables(physeq, errorIfNULL=TRUE)
#' 
#' @param physeq (Required). \code{\link{sampleData-class}}, or \code{\link{phyloseq-class}}.
#'
#' @param errorIfNULL (Optional). Logical. Should the accessor stop with 
#'  an error if the slot is empty (\code{NULL})? Default \code{TRUE}.
#'
#' @return Character vector. The names of the variables in the sampleData
#'  data.frame. Essentially the column names. Useful for selecting model 
#'  and graphics parameters that interact with sampleData.
#' 
#' @seealso getSpecies species.names sample.names getTaxa
#'
#' @export
#'
#' @examples
#' data(enterotype)
#' sample.variables(enterotype)
sample.variables <- function(physeq, errorIfNULL=TRUE){
	colnames(sampleData(physeq, errorIfNULL))
}
################################################################################
#' Get the values for a particular variable in sampleData
#'
#' This is a simple accessor function for streamlining access
#' to values/vectors/factors/etc contained in the sampleData.
#'
#' @usage getVariable(physeq, varName)
#' 
#' @param physeq (Required). \code{\link{sampleData-class}}, or \code{\link{phyloseq-class}}.
#'
#' @param varName (Required). Character string of the variable name in \code{sampleData}.
#'  Use \code{sample.variables(physeq)} for available variables in your object.
#'
#' @return Data. The clas of the data depends on what the contents of sampleData.
#' 
#' @seealso getSpecies species.names sample.names getTaxa 
#'  \code{\link{sample.variables}}
#'
#' @export
#'
#' @examples
#' # Load the GlobalPatterns dataset into the workspace environment
#' data(GlobalPatterns)
#' # Look at the different values for SampleType 
#' getVariable(GlobalPatterns, "SampleType")
getVariable <- function(physeq, varName){
	if( is.null(sampleData(physeq, FALSE)) ){
		stop("Your phyloseq data object does not have a sample-data component\n",
			"Try ?sampleData for more details.")
	}
	return( as(sampleData(physeq), "data.frame")[, varName] )
}
################################################################################