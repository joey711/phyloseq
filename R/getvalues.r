################################################################################
#' Returns the abundance values of species \code{i} for 
#' all samples in \code{x}.
#'
#' This is a simple accessor function for investigating 
#' a single species-of-interest. 
#'
#' @usage getSamples(x, i)
#' @param x otuTable, or H.O. object containing an otuTable.
#' @param i a single taxa/species/OTU ID for which you want
#'  to know the abundance in each sample.
#'
#' @return An integer vector of the abundance values for 
#' each sample in \code{x} for species \code{i}
#' 
#' @seealso getSpecies species.names sample.names
#'
#' @rdname getSamples-methods
#' @docType methods
#' @export
#'
#' @examples #
setGeneric("getSamples", function(x, i) standardGeneric("getSamples"))
################################################################################
#' @aliases getSamples,otuTable-method
#' @rdname getSamples-methods
setMethod("getSamples", "otuTable", function(x, i){
	if( speciesAreRows(x) ){
		as(x, "matrix")[i, ]
	} else {
		as(x, "matrix")[, i]
	}
})
################################################################################
#' @aliases getSamples,phyloseq-method
#' @rdname getSamples-methods
setMethod("getSamples", "phyloseq", function(x, i){
	getSamples(otuTable(x), i)
})
################################################################################
#' Returns the abundance values of sample \code{i} for 
#' all species in \code{x}.
#'
#' This is a simple accessor function for investigating 
#' a single sample-of-interest. 
#'
#' @usage getSpecies(x, i)
#' @param x otuTable, or H.O. object containing an otuTable.
#' @param i a single sample for which you want
#'  to know the abundance of each species.
#'
#' @return An integer vector of the abundance values for 
#' each species in \code{x} for sample \code{i}
#' 
#' @seealso getSpecies species.names sample.names
#'
#' @rdname getSpecies-methods
#' @docType methods
#' @export
#'
#' @examples #
setGeneric("getSpecies", function(x, i) standardGeneric("getSpecies"))
#' @aliases getSpecies,otuTable-method
#' @rdname getSpecies-methods
setMethod("getSpecies", "otuTable", function(x, i){
	if( speciesAreRows(x) ){
		as(x, "matrix")[, i]
	} else {
		as(x, "matrix")[i, ]
	}
})
#' @aliases getSpecies,phyloseq-method
#' @rdname getSpecies-methods
setMethod("getSpecies", "phyloseq", function(x, i){
	getSpecies(otuTable(x), i)
})
################################################################################