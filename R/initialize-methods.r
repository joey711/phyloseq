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
# Initialization methods for the higher-order classes.
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
setMethod("initialize", "otuTax", function(.Object, ...) {
	.Object <- callNextMethod()
	species <- intersect(rownames(.Object@taxTab), species.names(.Object))
	.Object@taxTab    <- taxTab(.Object)[species, ]
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