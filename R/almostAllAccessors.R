################################################################################
### Accessor / subset methods.
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
#' # OTU1 <- otuTable(matrix(sample(0:5,250,TRUE),25,10), speciesAreRows=TRUE)
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
setMethod("otuTable", "ANY", function(object, errorIfNULL=TRUE){
	access(object, "otuTable", errorIfNULL) 
})
# # # Instantiate an otuTable from a raw abundance matrix.
# # # 
# # # @param object An abundance table in the form of an integer matrix. 
# # #  \code{speciesAreRows} must be specified in \code{...}-argument .
# # # @param ... The additional named argument \code{speciesAreRows} must be 
# # #  provided. A logical.
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
#' Build or access sampleMap objects.
#'
#' \code{sampleMap()} is both a constructor and accessor method. When the
#' argument is a data.frame, sampleMap() will attempt to create and return an 
#' sampleMap-class object. In this case, the rows should be named to match the
#' \code{sample.names} of the other objects to which it will ultimately be paired.
#' Alternatively, if the argument is an object that 
#' contains a sampleMap, including a sampleMap-class object, then the 
#' corresponding sampleMap is returned, as-is.
#'
#' This is the main method suggested for constructing sampleMap objects from 
#' a data.frame of sample covariates. Each row represents a different sample.
#' It is also the suggested method for accessing/subsetting
#' the sampleMap from a more complex object.
#'
#' @usage sampleMap(object, errorIfNULL=TRUE)
#'
#' @param object An object among the set of classes defined by the phyloseq 
#' package that contain sampleMap.
#'
#' @param errorIfNULL (Optional). Logical. Should the accessor stop with 
#'  an error if the slot is empty (\code{NULL})? Default \code{TRUE}. 
#'
#' @return A sampleMap object. It is either grabbed from the relevant slot
#' if \code{object} is complex, or built anew if \code{object} is a data.frame
#' representing the sample covariates of an experiment.
#'
#' @seealso otuTable taxTab tre phyloseq sampleMap-class
#' @aliases sampleMap samplemap
#'
#' @rdname sampleMap-methods
#' @docType methods
#' @export
#'
#' @examples #
#' # OTU1 <- otuTable(matrix(sample(0:5,250,TRUE),25,10), speciesAreRows=TRUE)
#' # tax1 <- taxTab(matrix("abc", 30, 8))
#' # map1 <- data.frame( matrix(sample(0:3,250,TRUE),25,10), 
#' #   matrix(sample(c("a","b","c"),150,TRUE), 25, 6) ) 
#' # map1 <- sampleMap(map1)
#' # ex1 <- phyloseq(OTU1, map1, tax1)
#' # sampleMap(ex1)
setGeneric("sampleMap", function(object, errorIfNULL=TRUE) standardGeneric("sampleMap"))
#' @rdname sampleMap-methods
#' @aliases sampleMap,ANY-method
setMethod("sampleMap", "ANY", function(object){
	access(object, "sampleMap", errorIfNULL)
})
# constructor; for creating sampleMap from a data.frame
#' @rdname sampleMap-methods
#' @aliases sampleMap,data.frame-method
setMethod("sampleMap", "data.frame", function(object){
	# Make sure there are no phantom levels in categorical variables
	object <- reconcile_categories(object)
	# Want dummy samples index names if missing
	if( all(rownames(object) == as.character(1:nrow(object))) ){
		rownames(object) <- paste("sa", 1:nrow(object), sep="")
	}	
	new("sampleMap", object)
})
################################################################################
#' Access taxTab slot, or instantiate taxonomyTable-class.
#'
#' \code{taxTab()} is both a constructor and accessor method. When the
#' argument is a character matrix, taxTab() will attempt to create and return a 
#' taxonomyTable-class object. In this case, the rows should be named to match the
#' \code{species.names} of the other objects to which it will ultimately be paired.
#' Alternatively, if the argument is an object that 
#' contains a taxonomyTable, including a taxonomyTable-class object, then the 
#' corresponding taxonomyTable is returned, as-is.
#'
#' This is the main method suggested for constructing taxonomyTable objects from 
#' a character matrix of taxonomy classifications. Each row represents a 
#' different species.
#' It is also the suggested method for accessing/subsetting
#' the taxonomyTable from a more complex object. \code{taxTab} is the slot name
#' that holds the taxonomyTable-class object in a multi-component phyloseq
#' object.
#'
#' @usage taxTab(object, errorIfNULL=TRUE)
#'
#' @param object An object among the set of classes defined by the phyloseq 
#' package that contain taxonomyTable.
#'
#' @param errorIfNULL (Optional). Logical. Should the accessor stop with 
#'  an error if the slot is empty (\code{NULL})? Default \code{TRUE}.
#'
#' @return A taxonomyTable object. It is either grabbed from the relevant slot
#' if \code{object} is complex, or built anew if \code{object} is a 
#' character matrix representing the taxonomic classification of species in the
#' experiment.
#'
#' @seealso otuTable sampleMap tre phyloseq
#' @aliases taxTab taxtab
#'
#' @rdname taxTab-methods
#' @docType methods
#' @export
#'
#' @examples #
#' # tax1 <- taxTab(matrix("abc", 30, 8))
#' # taxTab(ex1)
#' # tax1
setGeneric("taxTab", function(object, errorIfNULL=TRUE) standardGeneric("taxTab"))
#' @rdname taxTab-methods
#' @aliases taxTab,ANY-method
setMethod("taxTab",  "ANY", function(object, errorIfNULL=TRUE){
	access(object, "taxTab", errorIfNULL)
})
# Constructor; for creating taxonomyTable from a matrix.
#' @rdname taxTab-methods
#' @aliases taxTab,matrix-method
setMethod("taxTab", "matrix", function(object){
	# Want dummy species/taxa index names if missing
	if(is.null(rownames(object))){
		rownames(object) <- paste("sp", 1:nrow(object), sep="")
	}
	if(is.null(colnames(object))){
		colnames(object) <- paste("ta", 1:ncol(object), sep="")
	}	
	new("taxonomyTable", object)
})
#' @rdname taxTab-methods
#' @aliases taxTab taxtab
#' @export
taxtab <- taxTab
################################################################################
#' Access tre slot, or check/coerce to phylo4 class.
#'
#' \code{tre()} is an accessor OR coercion method. This is the main method suggested 
#' for accessing/subsetting
#' the phylogenetic tree (phylo4 class) from a more complex object.
#' \code{tre} is the slot name for trees. 
#' that holds the phylo4-class object in a multi-component phyloseq-package
#' object.
#' 
#' Note that the tip labels should be named to match the
#' \code{species.names} of the other objects to which it is going to be paired.
#' The initialization methods for more complex objects automatically check for
#' exact agreement in the set of species described by the phlyogenetic tree 
#' and the other components (taxonomyTable, otuTable). 
#' They also trim accordingly. Thus, the tip.labels in a phylo object
#' must be named to match the
#' \code{species.names} of the other objects to which it will ultimately be paired.
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
#' @return A phylo4 object. It is grabbed from the tre-slot. If object
#'  is not a phylogenetic tree and does not
#'  contain a tre-slot, then NULL is returned. This is a convenience wrapper
#'  of the \code{\link{access}} function.
#'
#' @seealso otuTable sampleMap taxTab phyloseq
#' @export
#' @rdname tre-methods
#' @docType methods
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
#' @aliases speciesAreRows speciesarerows
#' @export
speciesarerows <- function(object){
	otuTable(object)@speciesAreRows
}
#' @aliases speciesAreRows speciesarerows
#' @export
speciesAreRows <- speciesarerows
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
#' @aliases nsamples,ANY-method
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
#' @rdname sample.names-methods
#' @aliases sample.names,sampleMap-method
setMethod("sample.names", "sampleMap", function(x) x@sample.names)
#' @rdname sample.names-methods
#' @aliases sample.names,otuTable-method
setMethod("sample.names", "otuTable", function(x) x@sample.names)
#' @rdname sample.names-methods
#' @aliases sample.names,phyloseq-method
setMethod("sample.names", "phyloseq", function(x){
	x@otuTable@sample.names 
})
#' @aliases sample.names
#' @export
sampleNames <- sample.names
