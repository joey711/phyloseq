##############################################################################
#' Extension of the (base) transpose function for otuTable
#'
#' Only the otuTable of a more complex object is transposed.	
#'
#' @usage t(x)
#'
#' @param x An \code{otuTable}, or a higher-order class that contains an \code{otuTable}.
#'
#' @return The class of the object returned by \code{t} matches
#' the class of the argument, \code{x}. The \code{otuTable} is now
#' transposed, and \code{speciesAreRows} is toggled.
#'
#' @name t
#' @rdname transpose-methods
#' @docType methods
#' @export
#' @examples #
setGeneric("t")
#' @aliases t,otuTable-method
#' @rdname transpose-methods
setMethod("t", signature("otuTable"), function(x){
	#new("otuTable", t(x@.Data), speciesAreRows = (!speciesAreRows(x)))
	x <- otuTable( t(as(x, "matrix")), speciesAreRows=(!speciesAreRows(x)) )
	return(x)
})
##############################################################################
#' @aliases t,phyloseqFather-method
#' @rdname transpose-methods
setMethod("t", signature("phyloseqFather"), function(x){
	x@otuTable <- t( otuTable(x) )
	return(x)
})
##############################################################################


