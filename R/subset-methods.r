################################################################################
# subsetting functions
# Without these, the default coerces to the base object (e.g. matrix or data.frame)
################################################################################
#' Extract parts of otuTable
#'
#' @export
#' @aliases [,otuTable,ANY,ANY,ANY-method
#' @rdname extract-methods
setMethod("[", "otuTable", function(x,i,j,...){
	newx <- callNextMethod(x@.Data,i,j,drop=FALSE,...)
	new("otuTable", x@speciesAreRows, newx)
})
################################################################################
#' extract parts of sampleMap
#'
#' @export
#' @aliases [,sampleMap,ANY,ANY,ANY-method
#' @rdname extract-methods
setMethod("[", "sampleMap", function(x,i,j,...){
	new("sampleMap", callNextMethod(data.frame(x),i,j,drop=FALSE,...))
})
################################################################################
#' extract parts of taxonomyTable
#'
#' @export
#' @aliases [,taxonomyTable,ANY,ANY,ANY-method
#' @rdname extract-methods
setMethod("[", "taxonomyTable", function(x,i,j,...){
	new("taxonomyTable", callNextMethod(x@.Data,i,j,drop=FALSE,...))
})
################################################################################
################################################################################
