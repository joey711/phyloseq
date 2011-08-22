################################################################################
# subsetting functions
# Want subset of our special tables to still be those classes after subset
# Without these, the default coerces to the base object (e.g. matrix or data.frame)
################################################################################
#' extract parts of otuTable
#'
#' @name [
#' @aliases [,otuTable-method
#' @docType methods
#' @rdname extract-methods
#'
setMethod("[", "otuTable", function(x,i,j,...){
	newx <- callNextMethod(x@.Data,i,j,drop=FALSE,...)
	new("otuTable", x@speciesAreRows, newx)
})
################################################################################
#' extract parts of sampleMap
#'
#' @name [
#' @aliases [,sampleMap-method
#' @docType methods
#' @rdname extract-methods
#'
setMethod("[", "sampleMap", function(x,i,j,...){
	new("sampleMap", callNextMethod(data.frame(x),i,j,drop=FALSE,...))
})
################################################################################
#' extract parts of taxonomyTable
#'
#' @name [
#' @aliases [,taxonomyTable-method
#' @docType methods
#' @rdname extract-methods
#'
setMethod("[", "taxonomyTable", function(x,i,j,...){
	new("taxonomyTable", callNextMethod(x@.Data,i,j,drop=FALSE,...))
})
################################################################################
################################################################################
