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
#' @param object An object among the set of classes defined by the phyloseq 
#' package that contain an otuTable.
#'
#' @param ... (optional) ignored unless \code{object} is a matrix, in which case
#'  \code{speciesAreRows} is required as a named argument.
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
setGeneric("otuTable", function(object, ...) standardGeneric("otuTable"))
################################################################################
# # #' Return the otuTable object as-is.
# # #'
# # #' @param object An otuTable-class object. The abundance table. 
#' @name otuTable
#' @docType methods
#' @aliases otuTable,otuTable-method
#' @rdname otuTable-methods
setMethod("otuTable", "otuTable", function(object) access(object, "otuTable") )
################################################################################
# # #' Return the otuTable component from a H.O. object.
# # #'
# # #' @param object Any \emph{phyloseq}-package object that contains an otuTable.
#' @name otuTable
#' @docType methods
#' @aliases otuTable,phyloseqFather-method
#' @rdname otuTable-methods
setMethod("otuTable", "phyloseqFather", function(object) access(object, "otuTable") )
################################################################################
# # #' Instantiate an otuTable from a raw abundance matrix.
# # #' 
# # #' @param object An abundance table in the form of an integer matrix. 
# # #'  \code{speciesAreRows} must be specified in \code{...}-argument .
# # #' @param ... The additional named argument \code{speciesAreRows} must be 
# # #'  provided. A logical.
#' @name otuTable
#' @docType methods
#' @aliases otuTable,matrix-method
#' @rdname otuTable-methods
setMethod("otuTable", "matrix", function(object, ...){
	new("otuTable", object, ...)
})
################################################################################
# # #' Convert to matrix, then instantiate otuTable.
# # #' 
#' @aliases otuTable,data.frame-method
#' @rdname otuTable-methods
setMethod("otuTable", "data.frame", function(object, ...){
	otuTable(as(object, "matrix"), ...)
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
#' @param object An object among the set of classes defined by the phyloseq 
#' package that contain sampleMap.
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
setGeneric("sampleMap", function(object) standardGeneric("sampleMap"))
#' @rdname sampleMap-methods
#' @aliases sampleMap,ANY-method
setMethod("sampleMap", "ANY", function(object){
	access(object, "sampleMap")
})
# constructor; for creating sampleMap from a data.frame
#' @rdname sampleMap-methods
#' @aliases sampleMap,data.frame-method
setMethod("sampleMap", "data.frame", function(object){
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
#' @param object An object among the set of classes defined by the phyloseq 
#' package that contain taxonomyTable.
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
setGeneric("taxTab", function(object) standardGeneric("taxTab"))
#' @rdname taxTab-methods
#' @aliases taxTab,ANY-method
setMethod("taxTab",  "ANY", function(object){
	access(object, "taxTab")
})
# Constructor; for creating taxonomyTable from a matrix.
#' @rdname taxTab-methods
#' @aliases taxTab,matrix-method
setMethod("taxTab", "matrix", function(object){
	new("taxonomyTable", object)
})
#' @rdname taxTab-methods
#' @aliases taxTab taxtab
#' @export
taxtab <- taxTab
################################################################################
#' Access tre slot.
#'
#' \code{tre()} is an accessor method. This is the main method suggested 
#' for accessing/subsetting
#' the phylogenetic tree (phylo4 class) from a more complex object.
#' \code{tre} is the slot name for trees. 
#' that holds the phylo4-class object in a multi-component phyloseq-package
#' object. 
#' 
#' Note that the tip.labels should be named to match the
#' \code{species.names} of the other objects to which it is paired. The 
#' initialization methods for more complex objects automatically check for
#' exact agreement in the set of species described by the phlyogenetic tree 
#' and the other components (taxonomyTable, otuTable). 
#' They also trim accordingly. Thus, the tip.labels in a phylo object
#' must be named to match the
#' \code{species.names} of the other objects to which it will ultimately be paired.
#' 
#' @param object An object among the set of classes defined by the phyloseq 
#' package that contain a phylogenetic tree. If object already is a phylogenetic
#' tree (a component data class), then it is returned as-is.
#'
#' @return A phylo object. It is grabbed from the tre-slot. If object does not
#'  contain a tre-slot, then NULL is returned. This is a convenience wrapper
#'  of the \code{\link{access}} function.
#'
#' @seealso otuTable sampleMap taxTab phyloseq
#' @export
tre <- function(object){
	access(object, "tre")
}
################################################################################
#' Access speciesAreRows slot from otuTable objects.
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
#' # library("picante")
#' # data("phylocom")
#' # tree <- phylocom$phylo
#' # nspecies(tree)
setGeneric("nspecies", function(object) standardGeneric("nspecies"))
#' @rdname nspecies-methods
#' @aliases nspecies,ANY-method
setMethod("nspecies", "ANY", function(object) object@nspecies)
#' @rdname nspecies-methods
#' @aliases nspecies,phyloseqFather-method
setMethod("nspecies", "phyloseqFather", function(object) nspecies(otuTable(object)) )
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
#' # library("picante")
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
setMethod("species.names", "otuTable", function(object) object@species.names)
#' @rdname species.names-methods
#' @aliases species.names,taxonomyTable-method
setMethod("species.names", "taxonomyTable", function(object) rownames(object) )
#' @rdname species.names-methods
#' @aliases species.names,phyloseqFather-method
setMethod("species.names", "phyloseqFather", function(object){
	otuTable(object)@species.names
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
#' # library("picante")
#' # data("phylocom")
#' # tree <- phylocom$phylo
#' # OTU1 <- otuTable(phylocom$sample, speciesAreRows=FALSE)
#' # nsamples(OTU1)
#' # physeq1 <- phyloseq(OTU1, tree)
#' # nsamples(physeq1)
setGeneric("nsamples", function(object) standardGeneric("nsamples"))
#' @rdname nsamples-methods
#' @aliases nsamples,ANY-method
setMethod("nsamples", "ANY", function(object) object@nsamples)
#' @rdname nsamples-methods
#' @aliases nsamples,phyloseqFather-method
setMethod("nsamples", "phyloseqFather", function(object){
	otuTable(object)@nsamples
})
################################################################################
#' Get the sample names of the samples described by an object.
#' 
#' This method works on otuTable and sampleMap objects, as well as the more
#' complex objects that represent samples in an experiment.
#'
#' @param object An object among the set of classes defined by the phyloseq 
#' package that represent samples.
#'
#' @return A character vector of the names of the samples in \code{object}.
#'
#' @seealso species.names nsamples
#' @aliases sample.names sampleNames
#'
#' @rdname sample.names-methods
#' @docType methods
#' @export
#'
#' @examples #
#' # library("picante")
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
#' @aliases sample.names,phyloseqFather-method
setMethod("sample.names", "phyloseqFather", function(x){
	x@otuTable@sample.names 
})
#' @aliases sample.names
#' @export
sampleNames <- sample.names
################################################################################
#' Create a subset object with just the \code{otuTable} and \code{sampleMap} slots.
#' 
#' @param object An object among the set of classes defined by the phyloseq 
#' package that contain both a sampleMap and an otuTable.
#'
#' @return A phyloseq-package object of class \code{otuSam}.
#'
#' @aliases otuSam otusam
#' @export
#' 
#' @examples #
#' ## data(ex1)
#' ## class(ex1)
#' ## otuSam(ex1)
otuSam <- function(object){
	phyloseq(otuTable(object), sampleMap(object))
}
#' @aliases otusam otuSam
#' @export
otusam <- otuSam
################################################################################
#' Create a subset object with just the \code{otuTable} and \code{tre} slots.
#' 
#' @param object An object among the set of classes defined by the phyloseq 
#' package that contain both a phylogenetic tree and an otuTable.
#'
#' @return A phyloseq-package object of class \code{otuTree}.
#'
#' @aliases otuTree otutree
#' @export
#' 
#' @examples #
#' # library("picante")
#' # data("phylocom")
#' # tree <- phylocom$phylo
#' # OTU1 <- otuTable(phylocom$sample, speciesAreRows=FALSE)
#' # ex1 <- phyloseq(OTU1, tree)
#' # class(ex1)
#' # otuTree(ex1)
otuTree <- function(object){
	phyloseq(otuTable(object), tre(object))
}
#' @aliases otuTree otutree
#' @export
otutree <- otuTree
################################################################################
#' Create a subset object with just the \code{otuTable} and \code{taxTab} slots.
#' 
#' This is a convenience wrapper of \code{\link{phyloseq}} that constructs
#' an otuTax-class object from a (equal or) more complex object that contains
#' an otuTable and taxTab slots, effectively creating a subsetted object.
#'
#' This is an accessor method only. For building an otuTax object from its 
#' components, the \code{\link{phyloseq}} function is suggested.
#'
#' @param object An object among the set of classes defined by the phyloseq 
#' package that contain a taxonomyTable, and an otuTable.
#'
#' @return A phyloseq-package object of class \code{otuSamTax}.
#'
#' @aliases otuTax otutax
#' @export
#' 
#' @examples #
otuTax <- function(object){
	phyloseq(otuTable(object), taxTab(object))
}
#' @aliases otuTax otutax
#' @export
otutax <- otuTax
################################################################################
#' Subset just the otuSamTree portion of a H.O. object.
#' 
#' @param object An object among the set of classes defined by the phyloseq 
#' package that contain a phylogenetic tree, sampleMap, and an otuTable.
#'
#' @return A phyloseq-package object of class \code{otuSamTree}.
#'
#' @aliases otuSamTree otuSamtree otusamtree
#' @export
#' 
#' @examples #
otuSamTree <- function(object){
	phyloseq(otuTable(object), sampleMap(object), tre(object))
}
#' @aliases otuSamTree otuSamtree otusamtree
#' @export
otusamtree <- otuSamTree
################################################################################
#' Subset just the otuSamTax portion of a H.O. object.
#' 
#' @param object An object among the set of classes defined by the phyloseq 
#' package that contain a taxonomyTable, sampleMap, and an otuTable.
#'
#' @return A phyloseq-package object of class \code{otuSamTax}.
#'
#' @aliases otuSamTax otuSamtax otusamtax
#' @export
#' 
#' @examples #
otuSamTax <- function(object){
	phyloseq(otuTable(object), sampleMap(object), taxTab(object))
}
#' @aliases otuSamTax otuSamtax otusamtax
#' @export
otusamtax <- otuSamTax
################################################################################
#' Show the component objects classes and slot names.
#'
#' There are no arguments to this function. It returns a named character
#' when called, which can then be used for tests of component data types, etc.
#' 
#' @return a character vector of the component objects classes, where each 
#' element is named by the corresponding slot name in the higher-order 
#' phyloseq objects (objects containing more than 1 phyloseq data object).
#' @export
#' @examples #
#' #get.component.classes()
get.component.classes <- function(){
	# define classes vector
	component.classes <- c("otuTable", "sampleMap", "phylo4", "phylo", "taxonomyTable")
	# the names of component.classes needs to be the slot names to match getSlots / splat
	names(component.classes) <- c("otuTable", "sampleMap", "tre", "old-tre", "taxTab")	
	return(component.classes)
}
################################################################################
#' Convert phyloseq objects into a named list of the component type (class)
#'
#' @param x An object of a class defined by the phyloseq-package. Component
#'  data and complex classes are both acceptable.
#' 
#' @return A named list, where each element is a component object that was contained 
#' in the argument, \code{x}. Each element is named for the object class it contains
#' If \code{x} is already a component data object,
#' then a list of length (1) is returned, also named.
#'
#' @seealso merge_phyloseq
#' @export
#' @examples #
splat.phyloseq.objects <- function(x){
	component.classes <- get.component.classes()
	# Check if class of x is among the component classes (not H.O.)
	if( class(x) %in% component.classes ){
		splatx <- list(x)
		names(splatx) <- names(component.classes)[component.classes==class(x)]
	} else {
		slotsx <- getSlots(class(x))
		splatx <- lapply(slotsx, function(iclass, slotsx, x){
			do.call(names(slotsx)[slotsx==iclass], list(x))
		}, slotsx, x)
	}
	return(splatx)
}
################################################################################
#' return the slot names of phyloseq objects.
#'
#' Like getSlots, but returns the class name if argument is component data object.
#' 
#' @param x An object of a class defined by the phyloseq-package. Component
#'  data and complex classes are both acceptable. 
#' 
#' @return identical to getSlots. A named character vector of the slot classes
#' of a particular S4 class, where each element is named by the slot name it
#' represents. If \code{x} is a component data object,
#' then a vector of length (1) is returned, named according to its slot name in
#' the higher-order objects.
#' 
#' @seealso merge_phyloseq
#' @export
#' @examples #
getslots.phyloseq <- function(x){
	# Check if class of x is among the component classes (not H.O.)
	component.classes <- get.component.classes()	
	if( class(x) %in% component.classes ){
		slotsx        <- as.character(class(x))
		names(slotsx) <- names(component.classes)[component.classes==class(x)]
	} else {
		slotsx <- getSlots(class(x))
	}
	return(slotsx)
}
################################################################################
#' General slot accessor function for phyloseq-package.
#'
#' This function is used internally by many convenience accessors and in 
#' many functions/methods that need to access a partiular component data type.
#' If something is wrong, or the slot is missing, the expected behavior is that
#' this function will return NULL. Thus, the output can be tested by 
#' \code{\link{is.null}} as verification of the presence of a particular 
#' data component.
#'
#' @param object An object of a class defined or extended by the phyloseq 
#'  package.
#'
#' @param slot A character string indicating the slot (not data class) of the 
#'  component data type that is desired.
#'
#' @return Returns the component object specified by the argument \code{slot}. 
#'  Returns NULL if slot does not exist. Returns \code{object} as-is 
#'  if it is a component class that already matches the slot name.
#' @seealso merge_phyloseq
#' @export
#' @examples #
#' ## data(ex1)
#' ## access(ex1, "taxTab")
#' ## access(ex1, "tre")
#' ## access(otuTable(ex1), "otuTable")
#' ## # Should return NULL:
#' ## access(otuTable(ex1), "sampleMap")
#' ## access(otuTree(ex1), "sampleMap")
#' ## access(otuSam(ex1), "tre")
access <- function(object, slot){
	component.classes <- get.component.classes()
	# Check if class of x is among the component classes (not H.O.)
	if( class(object) %in% component.classes ){
		# if slot-name matches object, return object as-is.
		if( component.classes[slot] == class(object) ){
			return( object )
		} else {
			return( NULL )
		}
	} else if(!slot %in% slotNames(object) ){
		return(NULL)
	} else {
		return( eval(parse(text=paste("object@", slot, sep=""))) )
	}
}
################################################################################
