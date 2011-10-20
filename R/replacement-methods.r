################################################################################
#' Assign to otuTable an object/value.
#'
#' @usage otuTable(x) <- value
#' @param x (Required). The object within which you will replace something
#' @param value (Required). The value with which you will replace that thing in x.
#'
#' @export
#' @docType methods
#' @rdname assign-otuTable
#' @aliases assign-otuTable otuTable<-
#'
#' @examples #
setGeneric("otuTable<-", function(x, value) standardGeneric("otuTable<-"))
#' @rdname assign-otuTable
#' @aliases otuTable<-,otuTable,otuTable-method
setMethod("otuTable<-", c("otuTable", "otuTable"), function(x, value){
	return(value)
})
#' @rdname assign-otuTable
#' @aliases otuTable<-,otuSam,otuTable-method
setMethod("otuTable<-", c("otuSam", "otuTable"), function(x, value){
	new("otuSam", otuTable=value, sampleMap=x@sampleMap)
})
#' @rdname assign-otuTable
#' @aliases otuTable<-,otuTree,otuTable-method
setMethod("otuTable<-", c("otuTree", "otuTable"), function(x, value){
	new("otuTree", otuTable=value, tre=x@tre)
})
#' @rdname assign-otuTable
#' @aliases otuTable<-,otuSamTax,otuTable-method
setMethod("otuTable<-", c("otuSamTax","otuTable"), function(x, value){
	new("otuSamTax", otuTable=value, sampleMap=x@sampleMap, taxTab=x@taxTab)
})
#' @rdname assign-otuTable
#' @aliases otuTable<-,otuSamTree,otuTable-method
setMethod("otuTable<-", c("otuSamTree", "otuTable"), function(x, value){
	new("otuSamTree", otuTable=value, sampleMap=x@sampleMap, tre=x@tre)
})
#' @rdname assign-otuTable
#' @aliases otuTable<-,otuSamTaxTree,otuTable-method
setMethod("otuTable<-", c("otuSamTaxTree", "otuTable"), function(x, value){
	new("otuSamTaxTree", otuTable=value,
		sampleMap=x@sampleMap, taxTab=x@taxTab, tre=x@tre)
})
################################################################################
#' Manually change speciesAreRows through assignment.
#'
#' The speciesAreRows slot is a logical indicating the orientation of the
#' abundance table contained in object \code{x}.
#'
#' @usage speciesarerows(x) <- value
#'
#' @param x An otuTable-class or higher-order object that contains an otuTable.
#'
#' @param value A logical of length equal to 1. If \code{length(value) > 1}, 
#'  the additional elements will be ignored. Only the first element is assigned
#'  to the speciesAreRows slot.
#'
#' @export
#' @docType methods
#' @rdname assign-speciesarerows
#' @aliases assign-speciesarerows speciesarerows<-
#'
#' @examples #
setGeneric("speciesarerows<-", function(x, value){
	standardGeneric("speciesarerows<-")
})
#' @rdname assign-speciesarerows
#' @aliases speciesarerows<-,otuTable,logical-method
setMethod("speciesarerows<-", c("otuTable", "logical"), function(x, value){
	x@speciesAreRows <- value[1]
	return(x)
})
#' @rdname assign-speciesarerows
#' @aliases speciesarerows<-,phyloseqFather,logical-method
setMethod("speciesarerows<-", c("phyloseqFather", "logical"), function(x, value){
	speciesAreRows(otuTable(x)) <- value
	return(x)
})
################################################################################
#' Assign species names to an object.
#'
#' A convenience / syntax wrapper around \code{\link{prune_samples}}.
#' This is typically used for subsetting an object to just those species
#' indicated by the character vector \code{value}.
#'
#' @usage species.names(x) <- value
#' @param x (Required). The object within which you will replace something
#' @param value (Required). The value with which you will replace that thing in x.
#'
#' @export
#' @rdname assign-species.names
#' @aliases assign-species.names species.names<-
#' @seealso prune_species
#' # # # data(ex1)
#' ### subset ex1 to just 5 random species
#' # # # species.names(ex1) <- sample(species.names(ex1), 5)
#' ### remove those species with only 1 or fewer individuals in ex1.
#' # # # species.names(ex1) <- species.names(ex1)[speciesSums(ex1) > 2]
"species.names<-" <- function(x, value){
	species <- intersect( value, species.names(x) )
	prune_species(species, x)
}
################################################################################
#' Assign sample names to an object.
#'
#' A convenience / syntax wrapper around \code{\link{prune_samples}}.
#' This is typically used for subsetting an object to just those samples
#' indicated by the character vector \code{value}.
#'
#' @usage sample.names(x) <- value
#' @param x (Required). The object within which you will replace something
#' @param value (Required). The value with which you will replace that thing in x.
#'
#' @export
#' @rdname assign-sample.names
#' @aliases assign-sample.names sample.names<-
#' @seealso \code{\link{prune_samples}}
#' @examples 
#' ## data(ex1)
#' # # # ex3 <- ex1
#' # # # sample.names(ex3) <- sample(sample.names(ex1), 5)
#' # # # sample.names(ex3) <- sample.names(ex3)[sampleMap(ex3)[,"Gender"]=="A"]
"sample.names<-" <- function(x, value){
	samples <- intersect( value, sample.names(x) )
	prune_samples(samples, x)
}
################################################################################
#' Assign to sampleMap an object/value.
#'
#' @usage sampleMap(x) <- value
#' @param x (Required). The object within which you will replace something
#' @param value (Required). The value with which you will replace that thing in x.
#'
#' @export
#' @rdname assign-sampleMap
#' @aliases assign-sampleMap sampleMap<-
#' @examples #
"sampleMap<-" <- function(x, value){
	x@sampleMap <- value
	sample.names(x) <- rownames(value)
	return(x)
}
################################################################################
#' Assign to taxTab an object/value.
#'
#' @usage taxTab(x) <- value
#' @param x (Required). The object within which you will replace something
#' @param value (Required). The value with which you will replace that thing in x.
#'
#' @export
#' @rdname assign-taxTab
#' @aliases assign-taxTab taxTab<-
#' @examples #
"taxTab<-" <- function(x, value){
	x@taxTab <- value
	species.names(x) <- rownames(value)
	return(x)
}
################################################################################
#' Assign to tre an object/value.
#'
#' This will automatically convert "phylo"-class trees (ape package) to
#' "phylo4"-class trees (phylobase package).
#'
#' @usage tre(x) <- value
#' @param x (Required). The object within which you will replace something
#' @param value (Required). The value with which you will replace that thing in x.
#'
#' @export
#' @docType methods
#' @rdname assign-tre
#' @aliases assign-tre tre<-
#' @examples #
setGeneric("tre<-", function(x, value) standardGeneric("tre<-"))
#' @rdname assign-tre
#' @aliases tre<-,otuTree,phylo-method
setMethod("tre<-", c("otuTree", "phylo"), function(x, value){
	x@tre <- as(value, "phylo4")
	species.names(x) <- value$tip.label
	return(x)
})
#' @rdname assign-tre
#' @aliases tre<-,otuTree,phylo4-method
setMethod("tre<-", c("otuTree", "phylo4"), function(x, value){
	x@tre <- value
	species.names(x) <- tipLabels(value) # value$tip.label
	return(x)
})
################################################################################
