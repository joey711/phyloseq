# Initialize packages. Copied approach from "S4 in 15 pages" near end.
#
# Since the goal of this code is to setup all the classes (and methods?) right away,
# it is also setup with it's name to be hidden from the user, accomplished by
# starting the function name with a dot. 
#
# Quote:
# "A relatively standard mechanism is to place all class, generic and method 
# definitions inside of a function. In the example below the function is called
# .initFoo 
# The reason it has a name that starts with a dot is so that the user will 
# not see it and inadvertantly call it.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
#
# Supposed to put all your classes and methods here within this
# dummy secret declaration function that gets hidden from user.
#.initphyloseq <- function(where) {
############################################################################
# Setup otuTable and its related methods 
############################################################################
#' The S4 class that holds taxa-abundance information.
#'
#' Because orientation of these tables can vary by method, the orientation is
#' defined explicitly in the \code{speciesAreRows} slot (a logical).
#' The \code{otuTable} class inherits the \code{\link{matrix}} class to store
#' abundance values.
#' Various standard subset and assignment nomenclature has been extended to apply
#' to the \code{otuTable} class, including square-bracket, \code{\link{t}}, etc.
#'
#' at-slot: speciesAreRows A single logical specifying the orientation of the 
#' abundance table.
#'
#' at-slot: nspecies A single positive integer specifying the number of distinct 
#' species/taxa.
#' 
#' at-slot: nsamples A single positive integer specifying the number of distinct
#' samples. Both \code{nsamples} and \code{nspecies} should not be assigned directly by user, but
#' are assessed internally during initialization of \code{otuTable} objects.
#' 
#' at-slot: species.names A character vector of the names of the taxa / species,
#' so that they can be accessed, modified, subsetted, without the need for
#' the user to consider the orientation of the abundance table.
#' 
#' at-slot: sample.names A character vector of the names of the samples, so that they
#' can be accessed, modified, subsetted, without the need for the user to consider
#' the orientation of the abundance table. Both \code{sample.names} and
#' \code{species.names} can
#' be assigned, e.g. \code{species.names(object)<-}, however, this has the effect
#' of internally replacing the otuTable to be a subset that matches the new name 
#' vector, rather than replacing the index names in-place.
#' 
#' at-slot: .Data This slot is inherited from the \code{\link{matrix}} class.
#'
#' @keywords internal
setClass("otuTable",
	representation(speciesAreRows="logical",
		nspecies="integer",
		nsamples="integer",
		species.names="character",
		sample.names="character"),
	contains = "matrix"
)
############################################################################
# Create an otuTable object from an abundance matrix.
# 
# (first unnamed entry), and a logical specifying the orientation in a named
# entry for \code{speciesAreRows}. Other slots in 
# the otuTable object
# are determined automatically from dimensions and dim-names. If no dim-names
# are provided, a default set is created with "sa" prefix for samples, and "sp"
# prefix for species / taxa. 
#
# This is an internal method. The suggested approach for creating an
# otuTable object from an abundance matrix is to use the 
# \link{\code{otuTable}} method instead.
# 
# @param .Object numeric (integer) matrix containing abundance data.
#
# @param speciesAreRows A single-element logical specifying the orientation
#  of the abundance table.
# 
# @param ... Optional additional arguments. Ignored.
# @seealso otuTable
# @keywords internal
# @rdname initialize-methods
# @examples 
# otuTable(matrix(sample(0:5,300,TRUE),30,10), speciesAreRows=FALSE)
setMethod("initialize", "otuTable", function(.Object, speciesAreRows, ...) { 
	.Object <- callNextMethod()
	.Object@speciesAreRows <- speciesAreRows
	if(speciesAreRows){
		# Want dummy species/sample index names if missing
		if(is.null(rownames(.Object@.Data))){
			rownames(.Object@.Data) <- paste("sp",1:nrow(.Object@.Data),sep="")
		}
		if(is.null(colnames(.Object@.Data))){
			colnames(.Object@.Data) <- paste("sa",1:ncol(.Object@.Data),sep="")
		}			
		.Object@nspecies <- nrow(.Object@.Data)
		.Object@species.names <- rownames(.Object@.Data)			
		.Object@nsamples <- ncol(.Object@.Data)
		.Object@sample.names <- colnames(.Object@.Data)
	} else {
		if(is.null(rownames(.Object@.Data))){
			rownames(.Object@.Data) <- paste("sa",1:nrow(.Object@.Data),sep="")
		}
		if(is.null(colnames(.Object@.Data))){
			colnames(.Object@.Data) <- paste("sp",1:ncol(.Object@.Data),sep="")
		} 
		.Object@nspecies      <- ncol(.Object@.Data)
		.Object@species.names <- colnames(.Object@.Data)
		.Object@nsamples      <- nrow(.Object@.Data)
		.Object@sample.names  <- rownames(.Object@.Data)			
	}
	.Object
})
############################################################################
#' The S4 object in the phyloseq package that holds sample data as a data.frame.
#'
#' at-slot: nsamples A single positive integer specifying the number of distinct samples.
#' \code{nsamples} should not be assigned directly by user, but
#' is assessed internally during initialization of \code{sampleMap} objects.
#'
#' at-slot: .Data data-frame data, inherited from the data.frame class.
#' 
#' at-slot: row.names Also inherited from the data.frame class;
#'  it should contain the sample names.
#' 
#' at-slot: names Inherited from the data.frame class.
#' 
#' at-slot: nsamples A single positive integer, indicating the number of samples
#'  (rows) represented by an object. Consistent with otuTable objects. 
#'  should not be assigned directly by user, but
#'  is assessed internally during instantiation of \code{sampleMap} objects.
#' 
setClass("sampleMap", representation(nsamples="integer"), contains="data.frame")
############################################################################
# Create an sampleMap object from a required data.frame class.
# (first unnamed entry). Other slots in \code{\linkS4class{sampleMap}} object
# are determined automatically from dimensions and row.names. If no dim-names
# are provided, a default set is created with "sa" prefix for samples.
# @param .Object numeric (integer) matrix containing abundance data.
# @rdname initialize-methods
setMethod("initialize", "sampleMap", function(.Object, ...) {
	.Object <- callNextMethod()
	.Object@nsamples <- nrow(.Object)
	# Want dummy samples index names if missing
	if( all(rownames(.Object) == as.character(1:nrow(.Object))) ){
		rownames(.Object) <- paste("sa", 1:nrow(.Object), sep="")
	}	
	.Object
})	
############################################################################
#' An S4 class that holds taxonomic classification data as a character
#' matrix.
#'
#' slots are:
#' 
#' nspecies -- A single positive integer specifying the number of 
#' distinct species/taxa.
#' \code{nspecies} should not be assigned directly by user, but
#' is assessed internally during instantiation of \code{taxonomyTable}
#' objects.
#' 
#' .Data -- This slot is inherited from the \code{\link{matrix}} class.
setClass("taxonomyTable", representation(nspecies="integer"), contains = "matrix")
############################################################################
# Create a \code{taxonomyTable} object from a character matrix.
#
# If \code{.Object} has no dim-names, then taxonomic-level and species/taxa names
# are added as dummy default.
#
# @param .Object character matrix containing taxonomic classifiers, where each row
# is a different species / taxa in a high-throughput sequencing experiment.
#
# @param ... Additional arguments. Ignored. 
#
# @rdname initialize-methods
# @examples
# taxtab1 <- new("taxonomyTable", matrix("abc", 5, 5))
setMethod("initialize", "taxonomyTable", function(.Object, ...) {
	.Object <- callNextMethod()
	.Object@nspecies <- nrow(.Object@.Data)
	# Want dummy species/taxa index names if missing
	if(is.null(rownames(.Object@.Data))){
		rownames(.Object@.Data) <- paste("sp", 1:nrow(.Object@.Data), sep="")
	}
	if(is.null(colnames(.Object@.Data))){
		colnames(.Object@.Data) <- paste("ta", 1:ncol(.Object@.Data), sep="")
	}
	.Object
})
############################################################################
# Define the higher-order classes in the phyloseq package.
# 
# 
# @import phylobase ape
# @importClassesFrom phylobase phylo4 phylo
########################################
#' Define the phyloseq class.
#'
#' Contains otuTable and sampleMap slots.
setClass("phyloseq", representation(otuTable="otuTable", sampleMap="sampleMap"))
########################################
#' Define the otuTree class.
#'
#' Contains otuTable, phylo (tre) slots.
#'
#' @importClassesFrom phylobase phylo
setClass("otuTree",  representation(otuTable="otuTable", tre="phylo"))				
########################################
#' Define the phyloseqTree class.
#'
#' Contains otuTable, sampleMap, phylo (tre) slots. Inherits the phyloseq and
#' otuTree classes.
#'
#' @importClassesFrom phylobase phylo
setClass("phyloseqTree", contains=c("phyloseq", "otuTree"))
########################################
#' Define the phyloseqTax class.
#'
#' Contains otuTable, sampleMap, and taxonomyTable (taxTab) slots. Inherits the
#' phyloseq class.
setClass("phyloseqTax", representation(taxTab="taxonomyTable"), contains="phyloseq")	
########################################
#' Define the phyloseqTaxTree class.
#'
#' Contains all (current) component classes: otuTable, sampleMap,
#' taxonomyTable (taxTab), and phylo (tre). 
#' Inherits phyloseqTax and phyloseqTree.
#' 
#' @importClassesFrom phylobase phylo
setClass("phyloseqTaxTree", contains=c("phyloseqTax", "phyloseqTree"))	
########################################
#' Define the otuTree4 class.
#'
#' Identical to otuTree class, but uses \code{phylo4} class for tree slot
#' instead of \code{phylo}.
#'
#' @importClassesFrom phylobase phylo4
setClass("otuTree4",  representation(otuTable="otuTable", tre="phylo4"))

############################################################################
# The initialization methods for the higher-order classes.
setMethod("initialize", "otuTree4", function(.Object, ...) {
	.Object <- callNextMethod()
	species <- intersect(species.names(.Object), tipLabels(.Object@tre)) #.Object@tre$tip.label)
	.Object@tre <- prune_species(species, .Object@tre)
	species.names(.Object@otuTable) <- species
	.Object
})
setMethod("initialize", "otuTree", function(.Object, ...) {
	.Object <- callNextMethod()
	species <- intersect(species.names(.Object), .Object@tre$tip.label)
	.Object@tre <- prune_species(species, .Object@tre)
	species.names(.Object@otuTable) <- species
	.Object
})
setMethod("initialize", "phyloseq", function(.Object, ...) {
	.Object <- callNextMethod()
	samples <- intersect(rownames(.Object@sampleMap), sample.names(.Object))
	.Object@sampleMap <- .Object@sampleMap[samples, ]
	.Object 
})
setMethod("initialize", "phyloseqTax", function(.Object, ...) {
	.Object <- callNextMethod()
	species <- intersect(rownames(.Object@taxTab), species.names(.Object))
	.Object@taxTab    <- taxTab(.Object)[species, ]
	species.names(.Object@otuTable) <- species
	samples <- intersect(rownames(.Object@sampleMap), sample.names(.Object))		
	.Object@sampleMap <- .Object@sampleMap[samples, ]		
	sample.names(.Object@otuTable) <- samples
	.Object 
})
setMethod("initialize", "phyloseqTree", function(.Object, ...) {
	.Object <- callNextMethod()
	species <- intersect(species.names(.Object), .Object@tre$tip.label)
	.Object@tre <- prune_species(species, .Object@tre)
	species.names(.Object@otuTable) <- species
	samples <- intersect(rownames(.Object@sampleMap), sample.names(.Object))		
	.Object@sampleMap <- .Object@sampleMap[samples, ]
	sample.names(.Object@otuTable)  <- samples		
	.Object 
})
setMethod("initialize", "phyloseqTaxTree", function(.Object, ...) {
	.Object <- callNextMethod()
	species <- intersect(rownames(.Object@taxTab), species.names(.Object))
	species <- intersect(species, .Object@tre$tip.label)
	.Object@taxTab    <- taxTab(.Object)[species, ]
	.Object@tre <- prune_species(species, .Object@tre)
	species.names(.Object@otuTable) <- species
	samples <- intersect(rownames(.Object@sampleMap), sample.names(.Object))		
	.Object@sampleMap <- .Object@sampleMap[samples, ]
	sample.names(.Object@otuTable)  <- samples
	.Object 
})

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
#' corresponding otuTable is returned, as-is.
#'
#' This is the main method suggested for constructing otuTable objects from 
#' an abundance matrix. It is also the suggested method for accessing subsetting
#' the otuTable from a more complex object.
#'
#' @param object An object among the set of classes defined by the phyloseq 
#' package that contain an otuTable.
#'
#' @param ... (optional) ignored unless object is a matrix. speciesAreRows is 
#' specified here. Must be a named argument.
#'
#' @return An otuTable object. It is either grabbed from the relevant slot
#' if \code{object} is complex, or built anew if \code{object} is an integer
#' matrix representing the species-abundance table.
#'
#' @seealso sampleMap taxTab tre phyloseq
#' @keywords OTU otuTable abundance
#' @aliases otutable
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
setMethod("otuTable", "otuTable", function(object) object)
setMethod("otuTable", "phyloseq", function(object) object@otuTable)
setMethod("otuTable", "otuTree", function(object) object@otuTable)
setMethod("otuTable", "otuTree4", function(object) object@otuTable)
# The following is for creating otuTables from a raw abundance matrix.
setMethod("otuTable", "matrix", function(object, ...){
	new("otuTable", object, ...)
})
# The following is for handling data.frames, as they should be
# converted to matrices
setMethod("otuTable", "data.frame", function(object, ...){
	otuTable(as(object, "matrix"), ...)
})
otutable <- otuTable
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
#' @seealso otuTable taxTab tre phyloseq
#' @keywords sampleMap covariates
#' @aliases sampleMap samplemap
#' @export
#' @examples #
#' # OTU1 <- otuTable(matrix(sample(0:5,250,TRUE),25,10), speciesAreRows=TRUE)
#' # tax1 <- taxTab(matrix("abc", 30, 8))
#' # map1 <- data.frame( matrix(sample(0:3,250,TRUE),25,10), 
#' #   matrix(sample(c("a","b","c"),150,TRUE), 25, 6) ) 
#' # map1 <- sampleMap(map1)
#' # ex1 <- phyloseq(OTU1, map1, tax1)
#' # sampleMap(ex1)
setGeneric("sampleMap", function(object) standardGeneric("sampleMap"))
setMethod("sampleMap", "phyloseq", function(object) object@sampleMap)
# The following is for creating sampleMap from a data.frame
setMethod("sampleMap", "data.frame", function(object){
	new("sampleMap", object)
})
samplemap <- sampleMap
################################################################################
#' Build or access taxonomyTable-class objects.
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
#' @keywords OTU tre phylo
#' @aliases taxTab taxtab
#' @export
#' @examples #
#' # OTU1 <- otuTable(matrix(sample(0:5,250,TRUE),25,10), speciesAreRows=TRUE)
#' # tax1 <- taxTab(matrix("abc", 30, 8))
#' # map1 <- data.frame( matrix(sample(0:3,250,TRUE),25,10), 
#' #   matrix(sample(c("a","b","c"),150,TRUE), 25, 6) ) 
#' # map1 <- sampleMap(map1)
#' # ex1 <- phyloseq(OTU1, map1, tax1)
#' # taxTab(ex1)
#' # tax1
setGeneric("taxTab", function(object) standardGeneric("taxTab"))
setMethod("taxTab", "phyloseqTax", function(object) object@taxTab)
# The following is for creating sampleMap from a data.frame
setMethod("taxTab", "matrix", function(object){
	new("taxonomyTable", object)
})
taxtab <- taxTab
taxonomyTable <- taxTab
taxonomytable <- taxTab
################################################################################
#' Access phylo-class objects.
#'
#' \code{tre()} is an accessor method. This is the main method suggested 
#' for accessing/subsetting
#' the phylogenetic tree (phylo class) from a more complex object.
#' \code{tre} is the slot name
#' that holds the phylo-class object in a multi-component phyloseq-package
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
#' package that contain a phylogenetic tree.
#'
#' @return A phylo object. It is grabbed from the tre slot
#'
#' @seealso otuTable sampleMap taxTab phyloseq
#' @keywords OTU tre tree phylo
#' @export
setGeneric("tre", function(object) standardGeneric("tre"))
setMethod("tre", "phylo", function(object) object)
setMethod("tre", "otuTree", function(object) object@tre)
setMethod("tre", "otuTree4", function(object) object@tre)
################################################################################
#' Access speciesAreRows slot from otuTable objects.
#'
#' @param object An object among the set of classes defined by the phyloseq 
#' package that are or contain an otuTable.
#'
#' @return A logical indicating the orientation of the otuTable.
#'
#' @seealso otuTable
#' @keywords OTU
#' @aliases speciesAreRows speciesarerows
#' @export
setGeneric("speciesarerows", function(object) standardGeneric("speciesarerows"))
setMethod("speciesarerows", "otuTable", function(object) object@speciesAreRows)
setMethod("speciesarerows", "otuTree", function(object){
	speciesarerows(otuTable(object))
})
setMethod("speciesarerows", "phyloseq", function(object){
	speciesarerows(otuTable(object))
})
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
#' @keywords taxa species OTU richness
#' @export
#' @examples #
#' # library("picante")
#' # data("phylocom")
#' # tree <- phylocom$phylo
#' # nspecies(tree)
setGeneric("nspecies", function(object) standardGeneric("nspecies"))	
setMethod("nspecies", "otuTable", function(object) object@nspecies)
setMethod("nspecies", "phylo", function(object) length(object$tip.label) )
setMethod("nspecies", "taxonomyTable", function(object) object@nspecies)
setMethod("nspecies", "otuTree", function(object) nspecies(otuTable(object)) )
setMethod("nspecies", "phyloseq", function(object) nspecies(otuTable(object)) )
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
#' @keywords taxa species OTU richness
#' @export
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
setMethod("species.names", "otuTable", function(object) object@species.names)
setMethod("species.names", "phylo", function(object) object$tip.label )
setMethod("species.names", "otuTree", function(object) species.names(otuTable(object)) )
setMethod("species.names", "otuTree4", function(object) species.names(otuTable(object)) )
setMethod("species.names", "phyloseq", function(object) species.names(otuTable(object)) )
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
#' @keywords samples
#' @export
#' @examples #
#' # library("picante")
#' # data("phylocom")
#' # tree <- phylocom$phylo
#' # OTU1 <- otuTable(phylocom$sample, speciesAreRows=FALSE)
#' # nsamples(OTU1)
#' # physeq1 <- phyloseq(OTU1, tree)
#' # nsamples(physeq1)
setGeneric("nsamples", function(object) standardGeneric("nsamples"))	
setMethod("nsamples", "otuTable", function(object) object@nsamples)
setMethod("nsamples", "otuTree", function(object) nsamples(otuTable(object)) )
setMethod("nsamples", "otuTree4", function(object) nsamples(otuTable(object)) )
setMethod("nsamples", "phyloseq", function(object) nsamples(otuTable(object)) )
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
#' @export
#' @examples #
#' # library("picante")
#' # data("phylocom")
#' # tree <- phylocom$phylo
#' # OTU1 <- otuTable(phylocom$sample, speciesAreRows=FALSE)
#' # sample.names(OTU1)
#' # ex1 <- phyloseq(OTU1, tree)
#' # sample.names(ex1)
setGeneric("sample.names", function(x) standardGeneric("sample.names"))	
setMethod("sample.names", "otuTable", function(x) x@sample.names)
setMethod("sample.names", "otuTree", function(x) sample.names(otuTable(x)) )
setMethod("sample.names", "otuTree4", function(x) sample.names(otuTable(x)) )
setMethod("sample.names", "phyloseq", function(x) sample.names(otuTable(x)) )
# Report an error if argument class is "NULL"
setMethod("sample.names", "NULL", function(x){
	print("argument class is NULL. Please check.")
})
sampleNames <- sample.names
################################################################################
#' Subset just the "otuTree" portion of a H.O. object.
#' 
#' @param object An object among the set of classes defined by the phyloseq 
#' package that contain both a phylogenetic tree and an otuTable.
#'
#' @return A phyloseq-package object of class \code{otuTree}.
#'
#' @seealso phyloseq
#' @keywords subset
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
setGeneric("otuTree", function(object) standardGeneric("otuTree"))
setMethod("otuTree", "otuTree", function(object){
	new("otuTree", otuTable=otuTable(object), tre=tre(object))
})
setMethod("otuTree", "phyloseqTree", function(object){ callNextMethod(object) })
setMethod("otuTree", "phyloseqTaxTree", function(object){ callNextMethod(object) })	
otutree <- otuTree
################################################################################
#' Subset just the "phyloseqTree" portion of a H.O. object.
#' 
#' @param object An object among the set of classes defined by the phyloseq 
#' package that contain a phylogenetic tree, sampleMap, and an otuTable.
#'
#' @return A phyloseq-package object of class \code{phyloseqTree}.
#'
#' @seealso phyloseq
#' @keywords subset
#' @aliases phyloseqTree phyloseqtree
#' @export
#' 
#' @examples #
setGeneric("phyloseqTree", function(object) standardGeneric("phyloseqTree"))
setMethod("phyloseqTree", "phyloseqTree", function(object){
	new("phyloseqTree", otuTable=otuTable(object), 
		sampleMap=sampleMap(object), tre=tre(object))		
})
setMethod("phyloseqTree", "phyloseqTaxTree",
	function(object){ callNextMethod(object) })
phyloseqtree <- phyloseqTree
################################################################################
#' Subset just the "phyloseqTax" portion of a H.O. object.
#' 
#' @param object An object among the set of classes defined by the phyloseq 
#' package that contain a taxonomyTable, sampleMap, and an otuTable.
#'
#' @return A phyloseq-package object of class \code{phyloseqTax}.
#'
#' @seealso phyloseq
#' @keywords subset
#' @aliases phyloseqTax phyloseqtax
#' @export
#' 
#' @examples #
setGeneric("phyloseqTax", function(object,...) standardGeneric("phyloseqTax"))
setMethod("phyloseqTax", "phyloseqTax", function(object){
	new("phyloseqTax", otuTable=otuTable(object), 
		sampleMap=sampleMap(object), taxTab=taxTab(object))
})
setMethod("phyloseqTax", "phyloseqTaxTree",
	function(object){ callNextMethod(object) })		
phyloseqtax <- phyloseqTax
################################################################################


#}
# Take this last line out when you are done testing...
#.initphyloseq()
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
