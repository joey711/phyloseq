################################################################################
#' Merge a subset of the species in \code{x} into one species/taxa/OTU.
#'
#' \code{mergespecies} is a method that takes as input an otuTable 
#' (or higher object) and a vector of species that should be merged.
#' It is intended to be able to operate at a low-level such that 
#' related methods, such as tipglom and taxglom can both reliably
#' call mergespecies for their respective purposes.
#'
#' @usage mergespecies(x, eqspecies, archetype=1)
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
#' @seealso \code{\link{tipglom}}, \code{\link{taxglom}}, \code{\link{merge_phyloseq}}
#'
#' @export
#' @docType methods
#' @rdname mergespecies-methods
#' @examples #
#' # # data(phylocom)
#' # # tree <- phylocom$phylo
#' # # otu  <- otuTable(phylocom$sample, speciesAreRows=FALSE)
#' # # otutree0 <- phyloseq(otu, tree)
#' # # plot(otutree0)
#' # # otutree1 <- mergespecies(otutree0, tree$tip.label[1:8], 2)
#' # # plot(otutree1)
setGeneric("mergespecies", function(x, eqspecies, archetype=1) standardGeneric("mergespecies"))
###############################################################################
# # #' Merge a subset of the taxa in an otuTable.
#'
#' @aliases mergespecies,otuTable-method
#' @rdname mergespecies-methods
setMethod("mergespecies", "otuTable", function(x, eqspecies, archetype=1){
	if( length(eqspecies) < 2 ){ return(x) }

	if( class(eqspecies) != "character" ){
		eqspecies <- species.names(x)[eqspecies]
	}
	# Shrink newx table to just those species in eqspecies
	newx <- prune_species(eqspecies, x)
	
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
	
	removeIndex <- which( species.names(x) %in% eqspecies[-keepIndex] )
	x <- prune_species(species.names(x)[-removeIndex], x)	
	return(x)
})
###############################################################################
# # #' Merge a subset of the taxa in a phylo-class tree object.
# # #'
# # #' @exportMethod mergespecies
# # #'
#' @aliases mergespecies,phylo-method
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
# # #' Merge a subset of the taxa in a phylo4-class tree object.
# # #'
# # #' @exportMethod mergespecies
# # #'
#' @aliases mergespecies,phylo4-method
#' @rdname mergespecies-methods
#'
#' @import phylobase
setMethod("mergespecies", "phylo4", function(x, eqspecies, archetype=1){
	y <- as(x, "phylo")
	y <- mergespecies(y, eqspecies, archetype)
	x <- as(y, "phylo4")
	return(x)
})
### phlyobase approach. Issues with complicated trees.
# # setMethod("mergespecies", "phylo4", function(x, eqspecies, archetype=1){
	# # if( length(eqspecies) < 2 ){ return(x) }

	# # if( class(eqspecies) != "character" ){
		# # eqspecies <- tipLabels(x)[eqspecies]  # x$tip.label[eqspecies]
	# # }
	# # if( class(archetype) != "character" ){
		# # keepIndex <- archetype
	# # } else {
		# # keepIndex <- which(eqspecies==archetype)
	# # }
	# # removeIndex <- which( tipLabels(x) %in% eqspecies[-keepIndex] )
	# # x           <- subset(x, tips.exclude=removeIndex)# drop.tip(x, removeIndex)
	# # return(x)
# # })
################################################################################
# # #' Merge a subset of the taxa in a complex phyloseq object.
# # #'
# # #' @exportMethod mergespecies
# # #'
#' @aliases mergespecies,phyloseqFather-method
#' @rdname mergespecies-methods
setMethod("mergespecies", "phyloseqFather", function(x, eqspecies, archetype=1){
	comp_list   <- splat.phyloseq.objects(x)
	merged_list <- lapply(comp_list, mergespecies, eqspecies, archetype)
	# the element names can wreak havoc on do.call
	names(merged_list) <- NULL
	# Re-instantiate the combined object
	do.call("phyloseq", merged_list)
})
###############################################################################
#' @aliases mergespecies,sampleMap-method
#' @rdname mergespecies-methods
setMethod("mergespecies", "sampleMap", function(x, eqspecies, archetype=1){
	return(x)
})
###############################################################################
#' @aliases mergespecies,taxonomyTable-method
#' @rdname mergespecies-methods
setMethod("mergespecies", "taxonomyTable", function(x, eqspecies, archetype=1){
	if( length(eqspecies) < 2 ){ return(x) }

	if( class(eqspecies) != "character" ){
		eqspecies <- species.names(x)[eqspecies]
	}
	
	if( class(archetype) != "character" ){
		keepIndex <- archetype
	} else {
		keepIndex <- which(eqspecies==archetype)
	}
	
	removeIndex <- which( species.names(x) %in% eqspecies[-keepIndex] )
	x <- prune_species(species.names(x)[-removeIndex], x)	
	return(x)
})
################################################################################
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
# tree = tre(ex4)
# tree1 = mergespecies(tree, tree$tip.label[1:500], 2)
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