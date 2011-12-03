################################################################################
#' Build or access otuTable objects.
#'
#' \code{otuTable()} is both a constructor and accessor method. When the first
#' argument is a matrix, otuTable() will attempt to create and return an 
#' otuTable-class object,
#' which further depends on whether or not speciesAreRows is provided as an
#' additional argument. Alternatively, if the first argument is an object that 
#' contains an otuTable, including an otuTable-class object, then the 
#' corresponding otuTable is returned, as the component object by itself.
#' This is a convenience wrapper on the more general \code{\link{access}} function
#' specific for grabbing the otuTable of an object.
#' It should work on both otuTable component objects, and higher-order classes
#' that contain an otuTable slot.
#'
#' This is the main method suggested for constructing otuTable objects from 
#' an abundance matrix. It is also the suggested method for accessing subsetting
#' the otuTable from a more complex object.
#'
#' @usage otuTable(object, speciesAreRows, errorIfNULL=TRUE)
#'
#' @param object (Required). A phyloseq object.
#'
#' @param speciesAreRows (Conditionally optional). Ignored unless 
#'  \code{object} is a matrix, in which case \code{speciesAreRows} is required.
#'
#' @param errorIfNULL (Optional). Logical. Should the accessor stop with 
#'  an error if the slot is empty (\code{NULL})? Default \code{TRUE}. Ignored
#'  if \code{object} argument is a matrix (constructor invoked instead).
#'
#' @return An otuTable object. It is either grabbed from the relevant slot
#' if \code{object} is complex, or built anew if \code{object} is an integer
#' matrix representing the species-abundance table. If \code{object} is a
#' data.frame, then an attempt is made to coerce it to an integer matrix and
#' instantiate an otuTable object.
#'
#' @name otuTable
#' @seealso sampleMap taxTab tre phyloseq
#' @aliases otuTable otutable otuTable,-method otuTable,otuTable-method
#' @docType methods
#' @rdname otuTable-methods
#' @export
#' @examples #
#' # OTU1 <- otuTable(matrix(sample(0:5,250,TRUE),25,10), TRUE)
#' # tax1 <- taxTab(matrix("abc", 30, 8))
#' # map1 <- data.frame( matrix(sample(0:3,250,TRUE),25,10), 
#' #    matrix(sample(c("a","b","c"),150,TRUE), 25, 6) )
#' # map1 <- sampleMap(map1) 
#' # ex1 <- phyloseq(OTU1, map1, tax1)
#' # otuTable(ex1)
setGeneric("otuTable", function(object, speciesAreRows, errorIfNULL=TRUE){
	standardGeneric("otuTable")	
})
# Access the otuTable slot, or return an otuTable as-is.
#' @aliases otuTable,phyloseq-method
#' @rdname otuTable-methods
setMethod("otuTable", "phyloseq", function(object, errorIfNULL=TRUE){
	access(object, "otuTable", errorIfNULL) 
})
# return the otuTable as-is, via access() for consistency.
#' @aliases otuTable,otuTable-method
#' @rdname otuTable-methods
setMethod("otuTable", "otuTable", function(object, errorIfNULL=TRUE){
	access(object, "otuTable", errorIfNULL) 
})
# Instantiate an otuTable from a raw abundance matrix.
#' @aliases otuTable,matrix-method
#' @rdname otuTable-methods
setMethod("otuTable", "matrix", function(object, speciesAreRows){
	# Want dummy species/sample index names if missing
	if(speciesAreRows){
		if(is.null(rownames(object))){
			rownames(object) <- paste("sp", 1:nrow(object), sep="")
		}
		if(is.null(colnames(object))){
			colnames(object) <- paste("sa", 1:ncol(object), sep="")
		}
	} else {
		if(is.null(rownames(object))){
			rownames(object) <- paste("sa",1:nrow(object),sep="")
		}
		if(is.null(colnames(object))){
			colnames(object) <- paste("sp",1:ncol(object),sep="")
		}
	}
	new("otuTable", object, speciesAreRows=speciesAreRows)
})
# # # Convert to matrix, then instantiate otuTable.
#' @aliases otuTable,data.frame-method
#' @rdname otuTable-methods
setMethod("otuTable", "data.frame", function(object, speciesAreRows){
	otuTable(as(object, "matrix"), speciesAreRows)
})
################################################################################
#' Returns the total number of individuals observed from each species.
#' 
#' A convenience function equivalent to rowSums or colSums, but where
#' the orientation of the otuTable is automatically handled. Can take 
#' more complex phyloseq objects, not just otuTable. Result always derived
#' from the abundance values in the associated otuTable, not other phyloseq
#' tables
#'
#' @usage speciesSums(x)
#'
#' @param x Any phyloseq object that is or contains an otuTable.
#' 
#' @return A named integer vector with length equal to the number of species
#'  in the table, name indicated the taxa ID, and value equal to the sum of
#'  all individuals observed for each taxa in \code{x}.
#'
#' @seealso sampleSums rowSums colSums
#' @export
#' @examples #
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
#' the orientation of the otuTable is automatically handled. Can take 
#' more complex phyloseq objects, not just otuTable. Result always derived
#' from the abundance values in the associated otuTable, not other phyloseq
#' tables.
#'
#' @usage sampleSums(x)
#'
#' @param x Any phyloseq-package object that is or contains an otuTable.
#' 
#' @return the total number of individuals present in each sample. A named 
#'  integer vector with length equal to the number of samples
#'  in the table, name indicated the sample ID, and value equal to the sum of
#'  all individuals observed for each sample in \code{x}.
#'
#' @seealso speciesSums rowSums colSums sum
#' @export
#' @examples #
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
