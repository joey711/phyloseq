####################################################################################
#' Transform the abundance count data in an \code{otuTable}, sample-by-sample.
#' 
#' This higher-order function transforms the sample counts of a species
#' abundance matrix according to a user-provided function or list of functions.
#' The counts of each sample will be transformed individually. No sample-sample 
#' interaction/comparison is possible by this method. 
#'
#' @usage transformsamplecounts(x, flist)
#'
#' @param x (Required). the \code{otuTable} or higher-order object that 
#'  contains an \code{otutable}.
#'
#' @param flist (Required). A list of functions that will be applied
#'  to the abundance counts of each sample.
#' 
#' @return The transformed \code{otuTable} or higher-order object with its
#'  \code{otuTable} transformed. In general, trimming is not expected by this 
#'  method, so it is suggested that the user provide only functions that return
#'  a full-length vector.
#'
#' @docType methods
#' @aliases transformsamplecounts TransformSampleCounts transformSampleCounts
#' @rdname transformcounts
#' @export
#'
#' @examples #
#' ## data(ex1)
#' ## ex1r <- transformsamplecounts(ex1, rank)
transformsamplecounts <- function(x, flist){
	if( speciesarerows(x) ){
		newx <- apply(as(otuTable(x), "matrix"), 2, flist)
	} else {
		newx <- apply(as(otuTable(x), "matrix"), 1, flist)
	}
	otuTable(x) <- otuTable(newx, speciesAreRows=speciesarerows(x))
	return(x)
}
####################################################################################
# # #' @aliases transformsamplecounts TransformSampleCounts transformSampleCounts
#' @docType methods
#' @rdname transformcounts
#' @export
TransformSampleCounts <- transformsamplecounts
####################################################################################
# # #' @aliases transformsamplecounts TransformSampleCounts transformSampleCounts
#' @docType methods
#' @rdname transformcounts
#' @export
transformSampleCounts <- transformsamplecounts
####################################################################################
