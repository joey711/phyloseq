################################################################################
# colname2hex
################################################################################
#' Convert color name to hex value.
#'
#' For converting color name to color hex and providing an optional
#' opacity parameter (alpha)
#'
#' @param colname A character vector of \code{R} color names.
#' 
#' @param alpha The standard alpha parameter specifying opacity/transparency.
#'
#' @return A color hex value of length equal to length of \code{colname}.
#'
#' @keywords internal 
#' @seealso col2rgb rgb2hsv colors
#' @examples #
colname2hex <- function(colname, alpha=1){
	if( length(colname) == 1){
		do.call("hsv", c(as.list(rgb2hsv(col2rgb(colname))[, 1]), alpha=alpha))
	} else if( length(colname) >1 ){
		sapply(colname, colname2hex, alpha)
	}	
}
################################################################################