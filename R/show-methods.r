############################################################################	
#' Create nice looking dimnames.
#'
#' This is internal. Not especially useful for others.
#'
#' @param z the table to be summarized.
#' @keywords internal
dimsum = function(z){
	if( length(z) < 6 ){ 
		return(z)
	} else {
		paste(paste( head(z, min(c(2,length(z))) ),sep="",collapse=", "),
			" ... ",
			  paste( tail(z, min(c(2,length(z))) ),sep="",collapse=", "),
		sep="",collapse="")
	}
}
############################################################################
#' print a summary of data tables to screen.
#'
#' This is internal. Not especially useful for others.
#'
#' @param y the table to be printed.
#' @param typelabel a character string indicating the class you want printed
#' @param dim1label a character string indicating the rows dim name.
#' @param dim2label a character string indicating the cols dim name.
#' @keywords internal
tableprint = function(y, 
	typelabel="table", dim1label="rownames", dim2label="colnames"){
	### Print OTU table features.
	if( any(dim(y))){
		cat( paste(typelabel," [",dim(y)[1]," by ",dim(y)[2],"]:",sep=""),fill=TRUE)
		cat( paste(dim1label,": ", sep="") )
		cat( dimsum(rownames(y)), fill=TRUE ) 
		cat( paste(dim2label,": ",sep="") ) 
		cat( dimsum(colnames(y)), fill=TRUE ) 
		# print tiny matrix part
		if( nrow(y) < 6 ){
			print( y[, 1:min(c(4, ncol(y))), drop=FALSE] )
		} else {
			print( head(y[, 1:min(c(4, ncol(y))), drop=FALSE], min(c(3, nrow(y)))) )
			cat("... ", fill=TRUE)
		}
	}			
}
############################################################################
#' @aliases show,otuTable-method
#' @rdname show-methods
setMethod("show", "otuTable", function(object){
	y <- as(object, "matrix")
	if( object@speciesAreRows ){
		tableprint(y, "OTU Table", "Species", "Samples")			
	} else {
		tableprint(y, "OTU Table", "Samples", "Species")
	}			 
})
############################################################################
#' @aliases show,sampleData-method
#' @rdname show-methods
setMethod("show", "sampleData", function(object){
	tableprint(object, "Sample Data", "Samples", "Variables")
})
############################################################################
#' @aliases show,taxonomyTable-method
#' @rdname show-methods
setMethod("show", "taxonomyTable", function(object){
	tableprint(as(object, "matrix"), "Taxonomy Table", "Species", "Taxonomic Rank")		
})
############################################################################
#' method extensions to show for phyloseq classes.
#'
#' See the general documentation of \code{\link[methods]{show}} method for
#' expected behavior. 
#'
#' @seealso \code{\link[methods]{show}}
#' 
#' @export
#' @aliases show,phyloseq-method
#' @docType methods
#' @rdname show-methods
#' @examples
#' # data(ex1)
#' # print(ex1)
#' # ex1
setMethod("show", "phyloseq", function(object){
	cat( "phyloseq-class experiment-level object \n", fill=TRUE )
	slot_list <- splat.phyloseq.objects(object)
	lapply(slot_list, function(i){show(i);cat("\n")})
})
############################################################################