################################################################################
#' Agglomerate taxa of the same type.
#'
#' This method merges species that are of the same taxaonomic level,
#' using \code{mergespecies} similar in approach to \code{tipglom}, but using
#' the categorical taxonomic levels rather than patristic distances on the
#' phylogenetic tree. 
#'
#' @param object An \code{otuTable} or more complex object that
#'  contains an \code{otuTable}.
#'
#' @param tax Either a \code{taxonomyTable} object, or alternatively, a 
#'  character vector specifying the desired taxonomic group of each taxa in 
#'  \code{object}. If \code{tax} is a character vector, it must have length equal
#'  to the (original) number of taxa in \code{object}, and each element must be
#'  named according to the taxa ID (that is, the result of 
#'  \code{species.names(object)}). If \code{tax} is a character vector, than
#'  the \code{taxlevel} argument is ignored. If \code{object} already contains
#'  a \code{taxonomyTable} component in its \code{taxTab} slot, then 
#'  the \code{tax} argument is ignored. 
#'
#' @param taxlevel A single-element character specifying the taxonomic level
#'  (column name)
#'  in \code{tax}, the \code{taxonomyTable}, that you want to agglomerate over.
#'  The default value is \code{"Phylum"}. Note that this default may
#'  agglomerate too broadly for a given experiment, and the user is strongly
#'  encouraged to try different taxonomic levels.
#' 
#' @return A taxonomically-agglomerated object with class matching
#' the class of \code{object}.
#'
#' @rdname taxglom-methods
#' @docType methods
#' @export
#'
#' @examples
#' ##
setGeneric("taxglom", function(object, tax=NULL, taxlevel="Phylum") standardGeneric("taxglom"))
#' @rdname taxglom-methods
#' @aliases taxglom,otuTable,taxonomyTable-method
setMethod("taxglom", c("otuTable", "taxonomyTable"), function(object, tax=NULL, taxlevel="Phylum"){
	# vectorize the taxonomy table.
	tax <- as(tax, "matrix")[, taxlevel]
	taxglom.internal(object, tax)
})
#' @rdname taxglom-methods
#' @aliases taxglom,otuTable,character-method
setMethod("taxglom", c("otuTable", "character"), function(object, tax=NULL, taxlevel="Phylum"){
	taxglom.internal(object, tax)
})
#' @rdname taxglom-methods
#' @aliases taxglom,otuTax,ANY-method
setMethod("taxglom", "otuTax", function(object, tax=NULL, taxlevel="Phylum"){
	# vectorize the taxonomy table.
	tax <- as(taxTab(object), "matrix")[, taxlevel]
	taxglom.internal(object, tax)
})
################################################################################
#' taxglom core internal function.
#'
#' taxglom.internal makes all the glomming happen, and delegates the
#' object-handling issues to \code{mergespecies()}.
#'
#' @param object the object on which agglomeration is to take place.
#'
#' @param tax for this internal function, tax must be a character vector.
#' \code{tax} is a vector in the core internal function.
#' 
#' @return the agglomerated object. Class matches argument \code{object}.
#'
#' @keywords internal
taxglom.internal <- function(object, tax){
	# Remove NAs and useless
	tax <- tax[ !tax %in% c("", " ", "\t") ]
	tax <- tax[ !is.na(tax) ]	
	
	# Define the species cliques to loop through
	spCliques     <- tapply(names(tax), factor(tax), list)
	
	# successively merge taxa in object.
	for( i in names(spCliques)){
		# print(i)
		object <- mergespecies(object, eqspecies=spCliques[[i]])
	}
	return(object)
}
################################################################################
# test <- taxglom.internal(ex1, as(taxTab(ex1), "matrix")[, "Phylum"])
# testvec <- as(taxTab(ex1), "matrix")[, "Phylum", drop=TRUE]
# tapply(names(testvec), factor(testvec), length)