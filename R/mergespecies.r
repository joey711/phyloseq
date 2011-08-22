################################################################################
#' Merge a subset of the species in \code{x} into one species/taxa/OTU.
#'
#' \code{mergespecies} is a method that takes as input an otuTable 
#' (or higher object) and a vector of species that should be merged.
#' It is intended to be able to operate at a low-level such that 
#' related methods, such as tipglom and taxglom can both reliably
#' call mergespecies for their respective purposes.
#'
#' @param x The otuTable or phylo object, or a higher-order object that
#'   contains an otuTable or tree.
#' 
#' @param eqspecies The species names, or indices, that should be merged together.
#'  If \code{length(eqspecies) < 2}, then the object \code{x} will be returned
#'  safely unchanged. 
#' 
#' @param archetype The index of \code{eqspecies} indicating the species 
#'   that should be kept (default is 1) to represent the summed/merged group
#'   of species/taxa/OTUs. 
#'   If archetype is not an index or index-name in \code{eqspecies}, the
#'   first will be used, and the value in archetype will be used 
#'   as the index-name for the new species.
#'
#' @return The object, \code{x}, in its original class, but with the specified
#'   species merged into one entry in all relevant components.
#'
#' @seealso tipglom taxglom merge_phyloseq
#'
#' @export
#' @rdname mergespecies-methods
#' @examples #
#' ## data(ex1)
#' ## makenetwork(otuTable(ex1), TRUE)
setGeneric("mergespecies", function(x, eqspecies, archetype=1) standardGeneric("mergespecies"))
###############################################################################
#' Merge a subset of the taxa in an otuTable.
#'
#' @exportMethod mergespecies
#'
#' @aliases mergespecies,otuTable-method
#' @docType methods
#' @rdname mergespecies-methods
setMethod("mergespecies", "otuTable", function(x, eqspecies, archetype=1){
	if( length(eqspecies) < 2 ){ return(x) }

	newx <-  x
	if( class(eqspecies) != "character" ){
		eqspecies <- species.names(x)[eqspecies]
	}
	# Shrink newx table to just those species in eqspecies
	species.names(newx) <- eqspecies
	
	if( class(archetype) != "character" ){
		keepIndex = archetype
	} else {
		keepIndex = which(eqspecies==archetype)
	}
	
	if( speciesAreRows(x) ){
		x[eqspecies[keepIndex], ] <- sampleSums(newx)
	} else {
		x[, eqspecies[keepIndex]] <- sampleSums(newx)	
	}
	removeIndex = which( species.names(x) %in% eqspecies[-keepIndex] )
	species.names(x) <- species.names(x)[-removeIndex]	
	return(x)
})
###############################################################################
#' Merge a subset of the taxa in a phylo-class tree object.
#'
#' @exportMethod mergespecies
#'
#' @aliases mergespecies,phylo-method
#' @docType methods
#' @rdname mergespecies-methods
setMethod("mergespecies", "phylo", function(x, eqspecies, archetype=1){
	if( length(eqspecies) < 2 ){ return(x) }

	if( class(eqspecies) != "character" ){
		eqspecies <- x$tip.label[eqspecies]
	}
	if( class(archetype) != "character" ){
		keepIndex <- archetype
	} else {
		keepIndex <- which(eqspecies==archetype)
	}
	removeIndex <- which( x$tip.label %in% eqspecies[-keepIndex] )
	x           <- ape::drop.tip(x, removeIndex)
	return(x)
})
###############################################################################
#' Merge a subset of the taxa in a phylo4-class tree object.
#'
#' @exportMethod mergespecies
#'
#' @aliases mergespecies,phylo4-method
#' @docType methods
#' @rdname mergespecies-methods
#'
#' @import phylobase
#'
setMethod("mergespecies", "phylo4", function(x, eqspecies, archetype=1){
	if( length(eqspecies) < 2 ){ return(x) }

	if( class(eqspecies) != "character" ){
		eqspecies <- tipLabels(x)[eqspecies]  # x$tip.label[eqspecies]
	}
	if( class(archetype) != "character" ){
		keepIndex <- archetype
	} else {
		keepIndex <- which(eqspecies==archetype)
	}
	removeIndex <- which( tipLabels(x) %in% eqspecies[-keepIndex] )
	x           <- subset(x, tips.exclude=removeIndex)# drop.tip(x, removeIndex)
	return(x)
})
###############################################################################
#' Merge a subset of the taxa in an otuTree4 object.
#'
#' Note to developers: The class name otuTree4 is a placeholder.
#' plan would be to migrate all ape-derived "phylo" classes to
#' \code{phylo4} of phylobase package (which has a namespace). Lots more potential
#' with phylo4, phylo4d, etc. E.g. merge could store the "lost"
#' tips as data (phylo4d) associated with each tip. More effective
#' than making separate list object to keep track.
#'
#' @exportMethod mergespecies
#'
#' @aliases mergespecies,otuTree4-method
#' @docType methods
#' @rdname mergespecies-methods
setMethod("mergespecies", "otuTree4", function(x, eqspecies, archetype=1){
	if( length(eqspecies) < 2 ){ return(x) }

  #stopifnot( class(eqspecies) == "character" )
	if( class(eqspecies) != "character" ){
		stop("For reliable behavior, eqspecies must be species names (character) \n 
		when merging species in otuTree object.")
	}
	newOTU <- mergespecies(otuTable(x), eqspecies, archetype)
	newtre <- mergespecies(tre(x), eqspecies, archetype)
	phyloseq(newOTU, newtre)
})
###############################################################################
#' Merge a subset of the taxa in an otuTree object.
#'
#' @exportMethod mergespecies
#'
#' @aliases mergespecies,otuTree-method
#' @docType methods
#' @rdname mergespecies-methods
setMethod("mergespecies", "otuTree", function(x, eqspecies, archetype=1){
	if( length(eqspecies) < 2 ){ return(x) }

	#stopifnot( class(eqspecies) == "character" )
	if( class(eqspecies) != "character" ){
		stop("For reliable behavior, eqspecies must be species names (character) \n 
		when merging species in otuTree object.")
	}
	newOTU <- mergespecies(otuTable(x), eqspecies, archetype)
	newtre <- mergespecies(tre(x), eqspecies, archetype)
	phyloseq(newOTU, newtre)
})
###############################################################################
#' Merge a subset of the taxa in an phyloseq object.
#'
#' @exportMethod mergespecies
#'
#' @aliases mergespecies,phyloseq-method
#' @docType methods
#' @rdname mergespecies-methods
setMethod("mergespecies", "phyloseq", function(x, eqspecies, archetype=1){
	if( length(eqspecies) < 2 ){ return(x) }

	otuTable(x) <- mergespecies(otuTable(x), eqspecies, archetype)
	return(x)
})
###############################################################################
#' Merge a subset of the taxa in an phyloseqTree object.
#'
#' @exportMethod mergespecies
#'
#' @aliases mergespecies,phyloseqTree-method
#' @docType methods
#' @rdname mergespecies-methods
setMethod("mergespecies", "phyloseqTree", function(x, eqspecies, archetype=1){
	if( length(eqspecies) < 2 ){ return(x) }

	y <- mergespecies(otuTree(x), eqspecies, archetype)
	x <- phyloseq(otuTable(y), tre(y), sampleMap(x))
	return(x)
})
###############################################################################
#' Merge a subset of the taxa in an phyloseqTaxTree object.
#'
#' @exportMethod mergespecies
#'
#' @aliases mergespecies,phyloseqTaxTree-method
#' @docType methods
#' @rdname mergespecies-methods
setMethod("mergespecies", "phyloseqTaxTree", function(x, eqspecies, archetype=1){
	if( length(eqspecies) < 2 ){ return(x) }

	y <- mergespecies(otuTree(x), eqspecies, archetype)
	x <- phyloseq(otuTable(y), tre(y), taxTab(x), sampleMap(x))
	return(x)
})
# ################################################################################
# # Example of higher-order phyloseq object species merge
# mergespecies(ex4, species.names(ex4)[1:5])
# mergespecies(phyloseqTree(ex4), species.names(ex4)[1:5])
# mergespecies(phyloseq(ex4), 1:5)
# ################################################################################
# # Example of otuTree species merge
# otutree  = otuTree(ex4)
# otutree1 = mergespecies(otutree, tre(otutree)$tip.label[1:9300])
# plot(tre(otutree1))
# # Not run, species indices not equivalent between phylo and otuTable.
# # Must use names (character):
# otutree2 = mergespecies(otutree, 1:9300)
################################################################################
# # Examples of tree species merge:
# ex1
# tree = tre(ex4)
# tree1 = mergespecies(tree, tree$tip.label[1:500], 2)

# data(phylocom)
# tree = phylocom$phylo
# tree1 = mergespecies(tree, tree$tip.label[1:8], 2)
# tree2 = mergespecies(tree, 12:15, 1)
# tree3 = mergespecies(tree, 12:15, 2)
# # Not run, won't know what merged:
# #tree3 = mergespecies(tree, sample(1:32,15), 2)
# par(mfcol=c(2,2))
# plot(tree,  main="Tree 0")
# plot(tree1, main="Tree 1")
# plot(tree2, main="Tree 2")
# plot(tree3, main="Tree 3")
################################################################################
# # Example of otuTable species merge
# x4 = otuTable(matrix(sample(0:15,100,TRUE),40,10), speciesAreRows=TRUE)
# mergespecies(x4, c("sp1", "sp3", "sp10", "sp28"), "sp10")
# mergespecies(x4, c("sp1", "sp3", "sp10", "sp28", "sp35") )
# mergespecies(x4, c("sp1", "sp3", "sp10", "sp28", "sp35") )@.Data
# mergespecies(x4, 5:25)@.Data
# Not run:
# mergespecies(x4, c("sp1", "sp3", "sp10", "sp28", "sp35", "") )
################################################################################