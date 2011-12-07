################################################################################
### Accessor / subset methods.
################################################################################
#' Access tre slot, or check/coerce to phylo4 class.
#'
#' \code{tre()} is an accessor OR coercion method. This is the main method suggested 
#' for accessing
#' the phylogenetic tree (\code{\link{phylo4-class}}) from a \code{\link{phyloseq-class}}.
#' Like other accessors (see See Also, below), the default behavior of this method
#' is to stop with an
#' error if \code{object} is a \code{phyloseq-class} but does not 
#' contain a phylogenetic tree.  
#' 
#' Note that the tip labels should be named to match the
#' \code{species.names} of the other objects to which it is going to be paired.
#' The initialization methods for more complex objects automatically check for
#' exact agreement in the set of species described by the phlyogenetic tree 
#' and the other components (taxonomyTable, otuTable). 
#' They also trim accordingly. Thus, the tip.labels in a phylo object
#' must be named to match the results of
#' \code{\link{species.names}} of the other objects to which it will ultimately be paired.
#'
#' @usage tre(object, errorIfNULL=TRUE)
#' 
#' @param object (Required). An instance of phyloseq-class
#'  that contains a phylogenetic tree. If object is a phylogenetic
#'  tree (a component data class), then it is returned as-is.
#'
#' @param errorIfNULL (Optional). Logical. Should the accessor stop with 
#'  an error if the slot is empty (\code{NULL})? Default \code{TRUE}.
#'
#' @return The \code{\link{phylo4-class}} object contained within \code{object};
#'  or NULL if \code{object} does not have a tree.
#'  This method stops with an error in the latter NULL case be default, which
#'  can be over-ridden by changing the value of \code{errorIfNULL} to \code{FALSE}.
#'
#' @seealso \code{\link{otuTable}}, \code{\link{sampleMap}}, \code{\link{taxTab}}
#'  \code{\link{phyloseq}}, \code{\link{merge_phyloseq}}
#' 
#' @export
#' @rdname tre-methods
#' @docType methods
#'
#' @examples
#' # data(ex1)
#' # tre(ex1)
setGeneric("tre", function(object, errorIfNULL=TRUE) standardGeneric("tre"))
#' @rdname tre-methods
#' @aliases tre,ANY-method
setMethod("tre", "ANY", function(object, errorIfNULL=TRUE){
	access(object, "tre", errorIfNULL)
})
# Constructor; for coercing "phylo4" from a "phylo".
#' @rdname tre-methods
#' @aliases tre,phylo-method
setMethod("tre", "phylo", function(object){
	as(object, "phylo4")
})
################################################################################
#' Access speciesAreRows slot from otuTable objects.
#'
#' @usage speciesarerows(object)
#'
#' @param object An object among the set of classes defined by the phyloseq 
#' package that are or contain an otuTable.
#'
#' @return A logical indicating the orientation of the otuTable.
#'
#' @seealso \code{\link{otuTable}}
#' @rdname speciesAreRows-methods
#' @docType methods
#' @export
#' @aliases speciesAreRows speciesarerows
setGeneric("speciesAreRows", function(object) standardGeneric("speciesAreRows"))
#' @rdname speciesAreRows-methods
#' @aliases speciesAreRows,ANY-method
setMethod("speciesAreRows", "ANY", function(object){NULL})
#' @rdname speciesAreRows-methods
#' @aliases speciesAreRows,otuTable-method
setMethod("speciesAreRows", "otuTable", function(object){object@speciesAreRows})
#' @rdname speciesAreRows-methods
#' @aliases speciesAreRows,phyloseq-method
setMethod("speciesAreRows", "phyloseq", function(object){
	speciesAreRows(otuTable(object))
})
#' @aliases speciesAreRows speciesarerows
#' @export
speciesarerows <- speciesAreRows
################################################################################
#' Get the number of taxa/species in an object.
#' 
#' This method works on otuTable, taxonomyTable objects, or the more
#' complex objects that represent species, as well as phylogenetic trees ``phylo''.
#'
#' @usage nspecies(object)
#'
#' @param object An object among the set of classes defined by the phyloseq 
#' package that represent species.
#'
#' @return An integer indicating the number of taxa / species.
#'
#' @seealso species.names
#'
#' @rdname nspecies-methods
#' @docType methods
#' @export
#'
#' @examples #
#' # # From "picante" package
#' # data("phylocom")
#' # tree <- phylocom$phylo
#' # nspecies(tree)
setGeneric("nspecies", function(object) standardGeneric("nspecies"))
#' @rdname nspecies-methods
#' @aliases nspecies,otuTable-method
setMethod("nspecies", "otuTable", function(object){
	if( speciesAreRows(object) ){
		return( length(rownames(object)) )
	} else {
		return( length(colnames(object)) )
	}
})
#' @rdname nspecies-methods
#' @aliases nspecies,taxonomyTable-method
setMethod("nspecies", "taxonomyTable", function(object){ nrow(object) })
#' @rdname nspecies-methods
#' @aliases nspecies,phyloseq-method
setMethod("nspecies", "phyloseq", function(object) nspecies(otuTable(object)) )
#' @rdname nspecies-methods
#' @aliases nspecies,phylo-method
setMethod("nspecies", "phylo", function(object) length(object$tip.label) )
#' @rdname nspecies-methods
#' @aliases nspecies,phylo4-method
#' @import phylobase
setMethod("nspecies", "phylo4", function(object) length(tipLabels(object)) )
################################################################################
#' Get the species / taxa names from an object.
#'
#' This method works on otuTable, taxonomyTable objects, or the more
#' complex objects that represent species, as well as phylogenetic trees
#' ``phylo''.
#'
#' @usage species.names(object)
#'
#' @param object An object among the set of classes defined by the phyloseq 
#' package that represent species.
#'
#' @return A character vector of the names of the species in \code{object}.
#'
#' @seealso nspecies
#'
#' @rdname species.names-methods
#' @docType methods
#' @export
#'
#' @examples #
#' # # From "picante" package
#' # data("phylocom")
#' # tree <- phylocom$phylo
#' # OTU1 <- otuTable(phylocom$sample, speciesAreRows=FALSE)
#' # species.names(tree)
#' # species.names(OTU1)
#' # physeq1 <- phyloseq(OTU1, tree)
#' # species.names(physeq1)
setGeneric("species.names", function(object) standardGeneric("species.names"))	
#' @rdname species.names-methods
#' @aliases species.names,otuTable-method
setMethod("species.names", "otuTable", function(object){
	if( speciesAreRows(object) ){
		return( rownames(object) )
	} else {
		return( colnames(object) )
	}
})
#' @rdname species.names-methods
#' @aliases species.names,taxonomyTable-method
setMethod("species.names", "taxonomyTable", function(object) rownames(object) )
#' @rdname species.names-methods
#' @aliases species.names,phyloseq-method
setMethod("species.names", "phyloseq", function(object){
	species.names(otuTable(object))
})
#' @rdname species.names-methods
#' @aliases species.names,sampleMap-method
setMethod("species.names", "sampleMap", function(object) NULL )
#' @rdname species.names-methods
#' @aliases species.names,phylo-method
setMethod("species.names", "phylo", function(object) object$tip.label )
#' @rdname species.names-methods
#' @aliases species.names,phylo4-method
#' @import phylobase
setMethod("species.names", "phylo4", function(object) tipLabels(object) )
################################################################################
#' Get the number of samples described by an object.
#' 
#' This method works on otuTable and sampleMap objects, as well as the more
#' complex objects that represent samples in an experiment.
#'
#' @usage nsamples(object)
#'
#' @param object An object among the set of classes defined by the phyloseq 
#' package that represent samples.
#'
#' @return An integer indicating the total number of samples.
#'
#' @seealso nspecies sample.names
#'
#' @rdname nsamples-methods
#' @docType methods
#' @export
#'
#' @examples #
#' # # From "picante" package
#' # data("phylocom")
#' # tree <- phylocom$phylo
#' # OTU1 <- otuTable(phylocom$sample, speciesAreRows=FALSE)
#' # nsamples(OTU1)
#' # physeq1 <- phyloseq(OTU1, tree)
#' # nsamples(physeq1)
setGeneric("nsamples", function(object) standardGeneric("nsamples"))
#' @rdname nsamples-methods
#' @aliases nsamples,otuTable-method
setMethod("nsamples", "otuTable", function(object){
	if( speciesAreRows(object) ){
		return( length(colnames(object)) )
	} else {
		return( length(rownames(object)) )
	}	
})
#' @rdname nsamples-methods
#' @aliases nsamples,phyloseq-method
setMethod("nsamples", "phyloseq", function(object){
	nsamples(otuTable(object))
})
################################################################################
#' Get the sample names of the samples described by an object.
#' 
#' This method works on otuTable and sampleMap objects, as well as the more
#' complex objects that represent samples in an experiment.
#'
#' @usage sample.names(x)
#'
#' @param x (Required). An object among the set of classes defined by the phyloseq 
#' package that represent samples.
#'
#' @return A character vector. The names of the samples in \code{x}.
#'
#' @seealso species.names nsamples
#' @aliases sample.names sampleNames
#'
#' @rdname sample.names-methods
#' @docType methods
#' @export
#'
#' @examples #
#' # # From "picante" package
#' # data("phylocom")
#' # tree <- phylocom$phylo
#' # OTU1 <- otuTable(phylocom$sample, speciesAreRows=FALSE)
#' # sample.names(OTU1)
#' # ex1 <- phyloseq(OTU1, tree)
#' # sample.names(ex1)
setGeneric("sample.names", function(x) standardGeneric("sample.names"))
# Unless otherwise specified, this should return a value of NULL
# That way, objects that do not explicitly describe samples all
# behave in the same (returning NULL) way.
#' @rdname sample.names-methods
#' @aliases sample.names,ANY-method
setMethod("sample.names", "ANY", function(x){ return(NULL) })
#' @rdname sample.names-methods
#' @aliases sample.names,sampleMap-method
setMethod("sample.names", "sampleMap", function(x) rownames(x) )
#' @rdname sample.names-methods
#' @aliases sample.names,otuTable-method
setMethod("sample.names", "otuTable", function(x){
	if( speciesAreRows(x) ){
		return( colnames(x) )
	} else {
		return( rownames(x) )
	}
})
#' @rdname sample.names-methods
#' @aliases sample.names,phyloseq-method
setMethod("sample.names", "phyloseq", function(x){
	sample.names(otuTable(x))
})
#' @aliases sample.names
#' @export
sampleNames <- sample.names
