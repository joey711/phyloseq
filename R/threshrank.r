####################################################################################
#' Thresholded rank transformation.
#' 
#' The lowest \code{thresh} values in \code{x} all get the value 'thresh'.
#'
#' @usage threshrank(x, thresh, keep0s=FALSE, ...)
#'
#' @param x The numeric vector to transform
#' @param thresh A single numeric value giving the threshold.
#' @param keep0s A logical determining whether 0's in \code{x} should remain 
#'  a zero-value in the output. If FALSE, zeros are treated as any other value.
#' @param ... Further arguments passes to the \code{\link{rank}} function.
#' 
#' @return A ranked, (optionally) thresholded numeric vector with length equal to
#'  \code{x}. Default arguments to \code{rank} are used, unless provided as
#'  additional arguments. 
#'
#' @seealso rank threshrankfun
#' @export 
#' @examples #
#' ## threshrank(sample(0:10, 100, TRUE), 50, keep0s=TRUE)
threshrank <- function(x, thresh, keep0s=FALSE, ...){
	if( keep0s ){ index0 <- which(x == 0) }
	x <- rank(x, ...)
	thresh <- thresh[1]
	x[x<thresh] <- thresh
	if( keep0s ){ x[index0] <- 0 }
	return(x)
}
####################################################################################
#' A closure version of the \code{threshrank} function.
#'
#' Takes the same arguments as \code{\link{threshrank}}, except for \code{x}, 
#' because the output is a function rather than a rank-transformed numeric. 
#' This approach is useful for creating an input to a higher-order function,
#' like "filterfun", that require a single-argument function as input. 
#'
#' @usage threshrankfun(thresh, keep0s=FALSE, ...)
#' @param thresh A single numeric value giving the threshold.
#' @param keep0s A logical determining whether 0's in \code{x} should remain 
#'  a zero-value in the output. If FALSE, zeros are treated as any other value.
#' @param ... Further arguments passes to the \code{\link{rank}} function.
#' 
#' @return A single-argument function with the options to threshrank set. 
#' @seealso topk
#' @export
#' @examples #
threshrankfun <- function(thresh, keep0s=FALSE, ...){
	function(x){
		threshrank(x, thresh, keep0s=FALSE, ...)
	}
}
####################################################################################
# example
# threshrank(sample(0:10, 100, TRUE), 50, keep0s=TRUE)
# transformsamplecounts( otuTable(ex4), threshrankfun(500))
#transformsamplecounts( otuTable(ex1), threshrankfun(500))
#transformsamplecounts( ex1, threshrankfun(500))