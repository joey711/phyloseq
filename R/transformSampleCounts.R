####################################################################################
#' Transform the abundance count data in an \code{otuTable}, sample-by-sample.
#' 
#' This higher-order function transforms the sample counts of a species
#' abundance matrix according to a user-provided function or list of functions.
#' The counts of each sample will be transformed individually. No sample-sample 
#' interaction/comparison is possible by this method. 
#'
#' @param x the \code{otuTable} or higher-order object that 
#'  contains an \code{otutable}.
#'
#' @param flist list of functions that will be applied to the abundance counts of
#'  each sample.
#' 
#' @return The transformed \code{otuTable} or higher-order object with its
#'  \code{otuTable} transformed. In general, trimming is not expected by this 
#'  method, so it is suggested that the user provide only functions that return
#'  a full-length vector.
#'
#' @seealso otuTable
#'
#' @export
#' @examples #
#' #data(ex1)
#' #ex1r <- transformsamplecounts(otuTable(ex1), rank)
setGeneric("transformsamplecounts", function(x, flist){
	standardGeneric("transformsamplecounts")})
setMethod("transformsamplecounts", "otuTable", function(x, flist){
	if( speciesarerows(x) ){
		newx <- apply(as(x, "matrix"), 2, flist)
	} else {
		newx <- apply(as(x, "matrix"), 1, flist)
	}
	new("otuTable", speciesarerows(x), newx)
})
setMethod("transformsamplecounts","otuTree",function(x, flist){
	otuTable(x) <- transformsamplecounts(otuTable(x), flist)
	return(x)
})
setMethod("transformsamplecounts","phyloseq",function(x, flist){
	otuTable(x) <- transformsamplecounts(otuTable(x), flist)
	return(x)
})
TransformSampleCounts <- transformsamplecounts
transformSampleCounts <- transformsamplecounts
####################################################################################
####################################################################################
