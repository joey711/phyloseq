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
# \code{\link{otuTable}} method instead.
# 
# @param .Object numeric (integer) matrix containing abundance data.
# 
# @param speciesAreRows A single-element logical specifying the orientation
#  of the abundance table.
#  
# @param ... Optional additional arguments. Ignored.
# 
# @seealso \code{\link{otuTable}}
#
# @examples 
# ## otuTable(matrix(sample(0:5,300,TRUE),30,10), speciesAreRows=FALSE)
#
#' @aliases initialize,otuTable-method
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
# Create an sampleMap object from a required data.frame class.
# (first unnamed entry). Other slots in \code{\link{sampleMap}} object
# are determined automatically from dimensions and row.names. If no dim-names
# are provided, a default set is created with "sa" prefix for samples.
# @param .Object numeric (integer) matrix containing abundance data.
# @seealso \code{\link{sampleMap}}
#
#' @aliases initialize,sampleMap-method
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
# @seealso \code{\link{taxTab}}
#
# @examples
# ## taxtab1 <- new("taxonomyTable", matrix("abc", 5, 5))
#
#' @aliases initialize,taxonomyTable-method
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
# Initialization methods for the higher-order classes.
############################################################################
# Need to add an automagic conversion for phylo -> phylo4
# Should be wrapped into the tre() build method, which would also be
# consistent with the other component class builders
#' @aliases initialize,otuTree-method
setMethod("initialize", "otuTree", function(.Object, ...) {
	.Object <- callNextMethod()
	.Object <- reconcile_species(.Object)
	.Object
})
#' @aliases initialize,otuTax-method
setMethod("initialize", "otuTax", function(.Object, ...) {
	.Object <- callNextMethod()
	.Object <- reconcile_species(.Object)
	.Object 
})
#' @aliases initialize,otuSam-method
setMethod("initialize", "otuSam", function(.Object, ...) {
	.Object <- callNextMethod()
	.Object <- reconcile_samples(.Object)
	.Object 
})
#' @aliases initialize,otuSamTax-method
setMethod("initialize", "otuSamTax", function(.Object, ...) {
	.Object <- callNextMethod()
	.Object <- reconcile_species(.Object)
	.Object <- reconcile_samples(.Object)
	.Object 
})
#' @aliases initialize,otuSamTree-method
setMethod("initialize", "otuSamTree", function(.Object, ...) {
	.Object <- callNextMethod()
	.Object <- reconcile_species(.Object)
	.Object <- reconcile_samples(.Object)
	.Object 
})
#' @aliases initialize,otuSamTaxTree-method
setMethod("initialize", "otuSamTaxTree", function(.Object, ...) {
	.Object <- callNextMethod()
	.Object <- reconcile_species(.Object)
	.Object <- reconcile_samples(.Object)
	.Object 
})
############################################################################
