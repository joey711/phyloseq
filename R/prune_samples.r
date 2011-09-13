################################################################################
#' Prune unwanted samples from a phyloseq object.
#' 
#' An S4 Generic method for removing (pruning) unwanted samples.
#'
#' @param samples A character vector of the samples in object x that you want to
#' keep. 
#'
#' @param x A phyloseq object.
#'
#' @return The class of the object returned by \code{prune_samples} matches
#' the class of the phyloseq object, \code{x}.
#'
#' @rdname prune_samples-methods
#' @docType methods
#' @export
#' @examples #
setGeneric("prune_samples", function(samples, x) standardGeneric("prune_samples"))
################################################################################
#' @aliases prune_samples,character,otuTable-method
#' @rdname prune_samples-methods
setMethod("prune_samples", signature("character", "otuTable"), function(samples, x){
	if( speciesarerows(x) ){
		x[, samples]
	} else {
		x[samples, ]
	}
})
################################################################################
#' @aliases prune_samples,character,sampleMap-method
#' @rdname prune_samples-methods
setMethod("prune_samples", signature("character", "sampleMap"), function(samples, x){
	x[samples, ]
})
################################################################################
#' @aliases prune_samples,character,otuSam-method
#' @rdname prune_samples-methods
setMethod("prune_samples", signature("character", "otuSam"), function(samples, x){
	x@sampleMap <- prune_samples(samples, x@sampleMap)
	x@otuTable  <- prune_samples(samples, x@otuTable)
	return(x)
})
################################################################################