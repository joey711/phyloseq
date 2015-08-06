################################################################################
# subsetting functions
# Without these, the default coerces to the base object (e.g. matrix or data.frame)
################################################################################
#' Method extensions to extraction operator for phyloseq objects.
#'
#' See the documentation for the \code{\link[base]{Extract}} generic,
#' defined in the R \code{\link[base]{base-package}}
#' for the expected behavior. 
#' 
#' One special exception to standard behavior of these methods in phyloseq is that
#' the \code{drop} argument is set internally to \code{FALSE}.
#' This helps avoid bugs during complicated subsetting with multiple components,
#' where it is necessary to be able to use a two dimensional indexing even
#' if one of those dimensions has only 1 rank. 
#' Put another way, these phyloseq-defined extractions never collapse their result
#' into a vector. See the documentation of \code{\link[base]{Extract}} for
#' more information about the \code{drop} argument.
#'
#' @param j See \code{\link[base]{Extract}}
#' 
#' @param ... See \code{\link[base]{Extract}}
#'
#' @seealso  \code{\link[base]{Extract}}
#' 
#' @export
#' 
#' @rdname extract-methods
#' @inheritParams base::Extract
#' @examples
#' data(esophagus)
#' nrow(otu_table(esophagus))
#' nrow(otu_table(esophagus)[1:5, ])
setMethod("[", "otu_table", function(x, i, j, ...){
	newx <- as(x, "matrix")[i, j, drop=FALSE]
	otu_table(newx, taxa_are_rows(x) )
})
# extract parts of sample_data
#
#' @export
#' @rdname extract-methods
setMethod("[", "sample_data", function(x, i, j, ...){
	sample_data( data.frame(x)[i, j, drop=FALSE] )
})
# extract parts of taxonomyTable
#
#' @export
#' @rdname extract-methods
setMethod("[", "taxonomyTable", function(x, i, j, ...){
  # Coerce to matrix, apply std extraction, reconstruct.
	return( tax_table(as(x, "matrix")[i, j, drop=FALSE]) )
})
# A numeric extraction method is already defined in Biostrings for XStringSet
# Add name-character-based extraction method for XStringSet
#
#' @importClassesFrom Biostrings XStringSet
#' @export
#' @rdname extract-methods
setMethod("[", c("XStringSet", "character"), function(x, i){
	index_vector = match(i, names(x), nomatch=NA_integer_)
	index_vector = index_vector[!is.na(index_vector)]
	if( length(index_vector) <= 0 ){
		warning("[,XStringSet: no valid seq-indices provided, NULL returned")
		return(NULL)
	}	
	if( length(index_vector) < length(i) ){
		warning("[,XStringSet: some seq-name indices invalid, omitted.")
	}
	# index_vector is an integer, subsetting now dispatches to standard
	x = x[index_vector]
	return(x)
})
################################################################################
################################################################################
