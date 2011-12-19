################################################################################
#' Build or access the otuTable.
#'
#' This is the suggested method for both constructing and accessing
#' Operational Taxonomic Unit (OTU) abundance (\code{\link{otuTable-class}}) objects.
#' When the first
#' argument is a matrix, otuTable() will attempt to create and return an 
#' otuTable-class object,
#' which further depends on whether or not \code{speciesAreRows} is provided as an
#' additional argument. 
#' Alternatively, if the first argument is an experiment-level (\code{\link{phyloseq-class}})
#' object, then the corresponding \code{otuTable} is returned.
#'
#' @usage otuTable(object, speciesAreRows, errorIfNULL=TRUE)
#'
#' @param object (Required). An integer matrix, \code{\link{otuTable-class}},
#'  or \code{\link{phyloseq-class}}.
#'
#' @param speciesAreRows (Conditionally optional). Logical; of length 1. Ignored
#'  unless \code{object} is a matrix, in which case it is is required.
#'
#' @param errorIfNULL (Optional). Logical. Should the accessor stop with 
#'  an error if the slot is empty (\code{NULL})? Default \code{TRUE}. Ignored
#'  if \code{object} argument is a matrix (constructor invoked instead).
#'
#' @return An \code{\link{otuTable-class}} object. 
#'
#' @seealso \code{\link{tre}}, \code{\link{sampleData}}, \code{\link{taxTab}}
#'  \code{\link{phyloseq}}, \code{\link{merge_phyloseq}}
#' 
#' @docType methods
#' @rdname otuTable-methods
#' @export
#' @examples #
#' # data(ex1)
#' # otuTable(ex1)
setGeneric("otuTable", function(object, speciesAreRows, errorIfNULL=TRUE){
	standardGeneric("otuTable")	
})
# Access the otuTable slot.
#' @aliases otuTable,phyloseq-method
#' @rdname otuTable-methods
setMethod("otuTable", "phyloseq", function(object, errorIfNULL=TRUE){
	access(object, "otuTable", errorIfNULL) 
})
# return the otuTable as-is.
#' @aliases otuTable,otuTable-method
#' @rdname otuTable-methods
setMethod("otuTable", "otuTable", function(object, errorIfNULL=TRUE){ return(object) })
# Instantiate an otuTable from a raw abundance matrix.
#' @aliases otuTable,matrix-method
#' @rdname otuTable-methods
setMethod("otuTable", "matrix", function(object, speciesAreRows){
	# instantiate first to check validity
	otutab <- new("otuTable", object, speciesAreRows=speciesAreRows)
	# Want dummy species/sample index names if missing
	if(speciesAreRows){
		if(is.null(rownames(otutab))){
			rownames(otutab) <- paste("sp", 1:nrow(otutab), sep="")
		}
		if(is.null(colnames(otutab))){
			colnames(otutab) <- paste("sa", 1:ncol(otutab), sep="")
		}
	} else {
		if(is.null(rownames(otutab))){
			rownames(otutab) <- paste("sa",1:nrow(otutab),sep="")
		}
		if(is.null(colnames(otutab))){
			colnames(otutab) <- paste("sp",1:ncol(otutab),sep="")
		}
	}
	return(otutab)
})
# # # Convert to matrix, then dispatch.
#' @aliases otuTable,data.frame-method
#' @rdname otuTable-methods
setMethod("otuTable", "data.frame", function(object, speciesAreRows){
	otuTable(as(object, "matrix"), speciesAreRows)
})
################################################################################
#' Returns the total number of individuals observed from each species/taxa/OTU.
#' 
#' A convenience function equivalent to rowSums or colSums, but where
#' the orientation of the otuTable is automatically handled.
#'
#' @usage speciesSums(x)
#'
#' @param x \code{\link{otuTable-class}}, or \code{\link{phyloseq-class}}.
#' 
#' @return A \code{\link{numeric-class}} with length equal to the number of species
#'  in the table, name indicated the taxa ID, and value equal to the sum of
#'  all individuals observed for each taxa in \code{x}.
#'
#' @seealso \code{\link{sampleSums}}, \code{\link{rowSums}}, \code{\link{colSums}}
#' @export
#' @examples
#' data(enterotype)
#' speciesSums(enterotype)
#' data(esophagus)
#' speciesSums(esophagus)
speciesSums <- function(x){
	x <- otuTable(x)
	if( speciesAreRows(x) ){
		rowSums(x)
	} else {
		colSums(x)
	}
}
speciessums <- speciesSums
################################################################################
#' Returns the total number of individuals observed from each sample.
#' 
#' A convenience function equivalent to rowSums or colSums, but where
#' the orientation of the otuTable is automatically handled.
#'
#' @usage sampleSums(x)
#'
#' @param x \code{\link{otuTable-class}}, or \code{\link{phyloseq-class}}.
#' 
#' @return A named \code{\link{numeric-class}}
#'  length equal to the number of samples
#'  in the \code{x}, name indicating the sample ID, and value equal to the sum of
#'  all individuals observed for each sample in \code{x}.
#'
#' @seealso \code{\link{speciesSums}}, \code{\link{rowSums}}, \code{\link{colSums}}
#' @export
#' @examples
#' data(enterotype)
#' sampleSums(enterotype)
#' data(esophagus)
#' sampleSums(esophagus)
sampleSums <- function(x){
	x <- otuTable(x)
	if( speciesAreRows(x) ){
		colSums(x)
	} else {
		rowSums(x)
	}
}
samplesums <- sampleSums
################################################################################