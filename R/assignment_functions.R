################################################################################
#' Assign to otuTable an object/value.
#'
#' @rdname assign-otuTable
#' @aliases assign-otuTable otuTable<-
#' @examples #
setGeneric("otuTable<-", function(x,value) standardGeneric("otuTable<-"))	
setMethod("otuTable<-", c("phyloseq","otuTable"), function(x,value){
	new("phyloseq", otuTable=value, sampleMap=x@sampleMap)
})
setMethod("otuTable<-", c("otuTree","otuTable"), function(x,value){
	new("otuTree", otuTable=value, tre=x@tre)
})
setMethod("otuTable<-", c("otuTree4", "otuTable"), function(x,value){
  new("otuTree4", otuTable=value, tre=x@tre)
})
setMethod("otuTable<-", c("phyloseqTax","otuTable"), function(x,value){
	new("phyloseqTax", otuTable=value, sampleMap=x@sampleMap, taxTab=x@taxTab)
})
setMethod("otuTable<-", c("phyloseqTree","otuTable"), function(x,value){
	new("phyloseqTree", otuTable=value, sampleMap=x@sampleMap, tre=x@tre)
})
setMethod("otuTable<-", c("phyloseqTaxTree","otuTable"), function(x,value){
	new("phyloseqTaxTree", otuTable=value,
		sampleMap=x@sampleMap, taxTab=x@taxTab, tre=x@tre)
})
################################################################################
#' Manually change speciesAreRows through assignment.
#'
#' @rdname assign-speciesarerows
#' @aliases assign-speciesarerows speciesarerows<-
#' @examples #
setGeneric("speciesarerows<-", function(x,value){
	standardGeneric("speciesarerows<-")
})
setMethod("speciesarerows<-", c("otuTable","logical"), function(x,value){
	x@speciesAreRows <- value
	return(x)
})
################################################################################
#' Assign species names to an object.
#'
#' This is typically used for subsetting an object to just those species
#' indicated by the character vector \code{value}.
#'
#' @rdname assign-species.names
#' @aliases assign-species.names species.names<-
#' @seealso prune_species
#' @examples #
setGeneric("species.names<-", function(x,value){
	standardGeneric("species.names<-")
})
setMethod("species.names<-", c("otuTable","character"), function(x,value){
	species <- intersect( value, species.names(x) )
	if( speciesarerows(x) ){
		x[species, ]
	} else {
		x[, species]
	}	
})
setMethod("species.names<-", c("otuTree","character"), function(x,value){
	species.names(otuTable(x)) <- value
	return(x)
})
setMethod("species.names<-", c("otuTree4","character"), function(x,value){
  species.names(otuTable(x)) <- value
	return(x)
})
setMethod("species.names<-", c("phyloseq","character"), function(x,value){
	species.names(otuTable(x)) <- value
	return(x)	
})
################################################################################
#' Assign sample names to an object.
#'
#' This is typically used for subsetting an object to just those samples
#' indicated by the character vector \code{value}.
#'
#' @rdname assign-sample.names
#' @aliases assign-sample.names sample.names<-
#' @seealso prune.samples
#' @examples #
setGeneric("sample.names<-", function(x,value){
	standardGeneric("sample.names<-")
})
setMethod("sample.names<-", c("otuTable","character"), function(x,value){
	samples <- intersect( value, sample.names(x) )
	if( speciesarerows(x) ){
		x[, samples]
	} else {
		x[samples, ]
	}
})
setMethod("sample.names<-", c("otuTree","character"), function(x,value){
	sample.names(otuTable(x))<-value
	return(x)
})
setMethod("sample.names<-", c("otuTree4","character"), function(x,value){
	sample.names(otuTable(x))<-value
	return(x)
})
setMethod("sample.names<-", c("phyloseq","character"), function(x,value){
	sample.names(otuTable(x))<-value
	return(x)	
})
################################################################################
#' Assign to sampleMap an object/value.
#'
#' @rdname assign-sampleMap
#' @aliases assign-sampleMap sampleMap<-
#' @examples #
setGeneric("sampleMap<-", function(x,value) standardGeneric("sampleMap<-"))	
setMethod("sampleMap<-", c("phyloseq","sampleMap"), function(x,value){
	x@sampleMap <- value
	sample.names(x) <- rownames(value)
	return(x)
})
################################################################################
#' Assign to taxTab an object/value.
#'
#' @rdname assign-taxTab
#' @aliases assign-taxTab taxTab<-
#' @examples #
setGeneric("taxTab<-", function(x,value) standardGeneric("taxTab<-"))	
setMethod("taxTab<-", c("phyloseq","taxonomyTable"), function(x,value){
	x@taxTab <- value
	species.names(x) <- rownames(value)
	return(x)
})
################################################################################
#' Assign to tre an object/value.
#'
#' @rdname assign-tre
#' @aliases assign-tre tre<-
#' @examples #
setGeneric("tre<-", function(x,value) standardGeneric("tre<-"))	
setMethod("tre<-", c("otuTree","phylo"), function(x, value){
	x@tre <- value
	species.names(x) <- value$tip.label
	return(x)
})
setMethod("tre<-", c("otuTree4","phylo4"), function(x, value){
	x@tre <- value
	species.names(x) <- tipLabels(value) # value$tip.label
	return(x)
})
################################################################################
