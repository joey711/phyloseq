################################################################################
#' Build or access the otu_table.
#'
#' This is the suggested method for both constructing and accessing
#' Operational Taxonomic Unit (OTU) abundance (\code{\link{otu_table-class}}) objects.
#' When the first
#' argument is a matrix, otu_table() will attempt to create and return an 
#' otu_table-class object,
#' which further depends on whether or not \code{taxa_are_rows} is provided as an
#' additional argument. 
#' Alternatively, if the first argument is an experiment-level (\code{\link{phyloseq-class}})
#' object, then the corresponding \code{otu_table} is returned.
#'
#' @usage otu_table(object, taxa_are_rows, errorIfNULL=TRUE)
#'
#' @param object (Required). An integer matrix, \code{\link{otu_table-class}},
#'  or \code{\link{phyloseq-class}}.
#'
#' @param taxa_are_rows (Conditionally optional). Logical; of length 1. Ignored
#'  unless \code{object} is a matrix, in which case it is is required.
#'
#' @param errorIfNULL (Optional). Logical. Should the accessor stop with 
#'  an error if the slot is empty (\code{NULL})? Default \code{TRUE}. Ignored
#'  if \code{object} argument is a matrix (constructor invoked instead).
#'
#' @return An \code{\link{otu_table-class}} object. 
#'
#' @seealso \code{\link{phy_tree}}, \code{\link{sample_data}}, \code{\link{tax_table}}
#'  \code{\link{phyloseq}}, \code{\link{merge_phyloseq}}
#' 
#' @docType methods
#' @rdname otu_table-methods
#' @export
#' @examples #
#' # data(GlobalPatterns)
#' # otu_table(GlobalPatterns)
setGeneric("otu_table", function(object, taxa_are_rows, errorIfNULL=TRUE){
	standardGeneric("otu_table")	
})
# Access the otu_table slot.
#' @aliases otu_table,phyloseq-method
#' @rdname otu_table-methods
setMethod("otu_table", "phyloseq", function(object, errorIfNULL=TRUE){
	access(object, "otu_table", errorIfNULL) 
})
# return the otu_table as-is.
#' @aliases otu_table,otu_table-method
#' @rdname otu_table-methods
setMethod("otu_table", "otu_table", function(object, errorIfNULL=TRUE){ return(object) })
# Instantiate an otu_table from a raw abundance matrix.
#' @aliases otu_table,matrix-method
#' @rdname otu_table-methods
setMethod("otu_table", "matrix", function(object, taxa_are_rows){
	# instantiate first to check validity
	otutab <- new("otu_table", object, taxa_are_rows=taxa_are_rows)
	# Want dummy species/sample index names if missing
	if(taxa_are_rows){
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
#' @aliases otu_table,data.frame-method
#' @rdname otu_table-methods
setMethod("otu_table", "data.frame", function(object, taxa_are_rows){
	otu_table(as(object, "matrix"), taxa_are_rows)
})
# Any less-specific class, not inherited by those above.
#' @aliases otu_table,ANY-method
#' @rdname otu_table-methods
setMethod("otu_table", "ANY", function(object, errorIfNULL=TRUE){
  access(object, "otu_table", errorIfNULL) 
})
################################################################################
#' Returns the total number of individuals observed from each species/taxa/OTU.
#' 
#' A convenience function equivalent to rowSums or colSums, but where
#' the orientation of the otu_table is automatically handled.
#'
#' @usage taxa_sums(x)
#'
#' @param x \code{\link{otu_table-class}}, or \code{\link{phyloseq-class}}.
#' 
#' @return A \code{\link{numeric-class}} with length equal to the number of species
#'  in the table, name indicated the taxa ID, and value equal to the sum of
#'  all individuals observed for each taxa in \code{x}.
#'
#' @seealso \code{\link{sample_sums}}, \code{\link{rowSums}}, \code{\link{colSums}}
#' @export
#' @examples
#' data(enterotype)
#' taxa_sums(enterotype)
#' data(esophagus)
#' taxa_sums(esophagus)
taxa_sums <- function(x){
	x <- otu_table(x)
	if( taxa_are_rows(x) ){
		rowSums(x)
	} else {
		colSums(x)
	}
}
################################################################################
#' Returns the total number of individuals observed from each sample.
#' 
#' A convenience function equivalent to rowSums or colSums, but where
#' the orientation of the otu_table is automatically handled.
#'
#' @usage sample_sums(x)
#'
#' @param x \code{\link{otu_table-class}}, or \code{\link{phyloseq-class}}.
#' 
#' @return A named \code{\link{numeric-class}}
#'  length equal to the number of samples
#'  in the \code{x}, name indicating the sample ID, and value equal to the sum of
#'  all individuals observed for each sample in \code{x}.
#'
#' @seealso \code{\link{taxa_sums}}, \code{\link{rowSums}}, \code{\link{colSums}}
#' @export
#' @examples
#' data(enterotype)
#' sample_sums(enterotype)
#' data(esophagus)
#' sample_sums(esophagus)
sample_sums <- function(x){
	x <- otu_table(x)
	if( taxa_are_rows(x) ){
		colSums(x)
	} else {
		rowSums(x)
	}
}
################################################################################
