################################################################################
#' Merge arguments into one phyloseq object.
#'
#' Takes a comma-separated list of phyloseq objects as arguments,
#' and returns the most-comprehensive single phyloseq object possible.
#'
#' Higher-order objects can be created if arguments are appropriate component data
#' types of different
#' classes, and this should mirror the behavior of the \code{\link{phyloseq}} method,
#' which is the suggested method if the goal is simply to create a higher-order
#' phyloseq object from different data types (1 of each class) describing the same experiment.
#' 
#' By contrast, this method is intended for situations in which one wants to combine
#' multiple higher-order objects, or multiple core component data objects (e.g. more than one
#' \code{otuTable}) that should be combined into one object.
#'
#' Merges are performed by first separating higher-order objects into
#' a list of their component objects; then, merging any component objects of the same class
#' into one object according to the behavior desribed in \code{\link{merge_phyloseq_pair}};
#' and finally, building back up a merged-object according to the constructor
#' behavior of the \code{\link{phyloseq}} method. If the arguments contain only a single
#' component type -- several otuTable objects, for example -- then a single merged object
#' of that component type is returned.
#' 
#' @usage merge_phyloseq(...)
#'
#' @param ... a comma-separated list of phyloseq objects. 
#'
#' @return Merges are performed by first separating higher-order objects into
#' a list of their component objects; then, merging any component objects of the same class
#' into one object according to the behavior desribed in \code{\link{merge_phyloseq_pair}};
#' and finally, re-building a merged-object according to the constructor
#' behavior of the \code{\link{phyloseq}} method. If the arguments contain only a single
#' component type -- several otuTable objects, for example -- then a single merged object
#' of the relevant component type is returned.
#'
#' Merges between 2 or more tree objects are ultimately done using 
#' \code{\link{consensus}} from the ape package.
#' This has the potential to limit somewhat the final data object, because trees
#' don't merge with other trees in the same granular manner as data tables, and
#' ultimately the species/taxa in higher-order phyloseq objects will be clipped to
#' what is contained in the tree. If this an issue, the tree component should
#' be ommitted from the argument list.
#'
#' @export
#'
#' @examples #
#' ## # Make a random complex object
#' ## OTU1 <- otuTable(matrix(sample(0:5,250,TRUE),25,10), speciesAreRows=TRUE)
#' ## tax1 <- taxTab(matrix("abc", 30, 8))
#' ## map1 <- data.frame( matrix(sample(0:3,250,TRUE),25,10), 
#' ##   matrix(sample(c("a","b","c"),150,TRUE), 25, 6) ) 
#' ## map1 <- sampleData(map1)
#' ## ex1 <- phyloseq(OTU1, map1, tax1)
#' ## x <- ex1
#' ## x <- phyloseq(ex1)
#' ## y <- taxTab(ex1)
#' ## merge_phyloseq(x, y)
#' ## merge_phyloseq(y, y, y, y)
#' ## # Try the simple example dataset, ex1
#' ## data(ex1)
#' ## merge_phyloseq( otutree(ex1), sampleData(ex1))
#' ## merge_phyloseq( otuSam(ex1), tre(ex1), taxTab(ex1))
merge_phyloseq <- function(...){
	arguments <- list(...)
	# create list of all components of all objects
	comp.list <- list()
	for( i in 1:length(arguments) ){
		comp.list <- c(comp.list, splat.phyloseq.objects(arguments[[i]]))
	}
	# loop through each component type. Note, list names redundant. will use this
	merged.list <- list()
	for( i in unique(names(comp.list)) ){ #i="taxTab"
		# check if length 1, if so, cat to merged.list.
		i.list <- comp.list[names(comp.list)==i]
		if( length(i.list) == 1 ){
			merged.list <- c(merged.list, i.list)
		} else {
			# else, loop through each identically-named objects.
			x1 <- i.list[[1]]
			for( j in 2:length(i.list)){
				x1 <- merge_phyloseq_pair(x1, i.list[[j]])
			}
			x1 <- list(x1)
			names(x1)   <- i
			merged.list <- c(merged.list, x1)			
		}
	}
	# Remove names to avoid any conflicts with phyloseq(), which does not need named-arguments
	names(merged.list) <- NULL

	# Use do.call for calling this variable-length, variable-content argument list.
	return( do.call(phyloseq, merged.list) )
}
################################################################################
#' Merge pair of phyloseq component data objects of the same class.
#'
#' Internal S4 methods to combine pairs of objects of classes specified in the
#' phyloseq package. These objects must be component data of the same type 
#' (class). This is mainly an internal method, provided to illustrate how
#' merging is performed by the more general \code{\link{merge_phyloseq}} function.
#'
#' The \code{\link{merge_phyloseq}} function is recommended in general.
#' 
#' Special note: trees are merged using \code{\link[ape]{consensus}}.
#'
#' @usage merge_phyloseq_pair(x, y) 
#'
#' @param x A character vector of the species in object x that you want to
#' keep -- OR alternatively -- a logical vector where the kept species are TRUE, and length
#' is equal to the number of species in object x. If \code{species} is a named
#' logical, the species retained is based on those names. Make sure they are
#' compatible with the \code{species.names} of the object you are modifying (\code{x}). 
#'
#' @param y Any \code{phyloseq} object.
#'
#' @return A single component data object that matches \code{x} and \code{y} 
#' arguments. The returned object will 
#' contain the union of the species and/or samples of each. If there is redundant
#' information between a pair of arguments of the same class, the values in \code{x} are
#' used by default. Abundance values are summed for \code{otuTable} objects 
#' for those elements that describe the same species and sample in \code{x}
#' and \code{y}. 
#'
#' @seealso \code{\link{merge_phyloseq}} \code{\link{merge_species}}
#'
#' @rdname merge_phyloseq_pair-methods
#' @docType methods
#' @export
#'
#' @examples #
#' ## # merge two simulated otuTable objects.
#' ## x  <- otuTable(matrix(sample(0:5,200,TRUE),20,10), speciesAreRows=TRUE)
#' ## y  <- otuTable(matrix(sample(0:5,300,TRUE),30,10), speciesAreRows=FALSE)
#' ## xy <- merge_phyloseq_pair(x, y)
#' ## yx <- merge_phyloseq_pair(y, x)
#' ## # merge two simulated taxTab objects
#' ## x <- taxTab(matrix("abc", 20, 6))
#' ## y <- taxTab(matrix("def", 30, 8))
#' ## xy <- merge_phyloseq_pair(x, y)
#' ## # merge two simulated sampleData objects
#' ## x <- data.frame( matrix(sample(0:3,250,TRUE),25,10), 
#' ##   matrix(sample(c("a","b","c"),150,TRUE),25,6) )
#' ## x <- sampleData(x)
#' ## y <- data.frame( matrix(sample(4:6,200,TRUE),20,10), 
#' ##   matrix(sample(c("d","e","f"),120,TRUE),20,8) )
#' ## y <- sampleData(y)
#' ## merge_phyloseq_pair(x, y)
#' ## data.frame(merge_phyloseq_pair(x, y))
#' ## data.frame(merge_phyloseq_pair(y, x))
setGeneric("merge_phyloseq_pair", function(x, y) standardGeneric("merge_phyloseq_pair"))
################################################################################
#' @aliases merge_phyloseq_pair,otuTable,otuTable-method
#' @rdname merge_phyloseq_pair-methods
setMethod("merge_phyloseq_pair", signature("otuTable", "otuTable"), function(x, y){
	specRrowsx   <- speciesAreRows(x)
	new.sp.names <- union(species.names(x), species.names(y))
	new.sa.names <- union(sample.names(x), sample.names(y))
	
	# Create the empty new matrix structure
	newx <- matrix(0, nrow=length(new.sp.names), ncol=length(new.sa.names),
		dimnames=list(new.sp.names, new.sa.names))

	# assign a standard speciesAreRows orientation to TRUE for x and y
	if( !speciesAreRows(x) ){ x <- t(x) }
	if( !speciesAreRows(y) ){ y <- t(y) }

	# "merge" by addition.
	newx[rownames(x), colnames(x)] <- x
	newx[rownames(y), colnames(y)] <- newx[rownames(y), colnames(y)] + y

	# Create the new otuTable object
	newx <- otuTable(newx, speciesAreRows=TRUE)

	# Return the orientation that was in x
	if( !specRrowsx ){ newx <- t(newx) }
	return(newx)
})
################################################################################
#' @aliases merge_phyloseq_pair,taxonomyTable,taxonomyTable-method
#' @rdname merge_phyloseq_pair-methods
setMethod("merge_phyloseq_pair", signature("taxonomyTable", "taxonomyTable"), function(x, y){
	new.sp.names <- union(rownames(x), rownames(y))
	new.ta.names <- union(colnames(x), colnames(y))
	
	# Create the empty new matrix structure
	newx <- matrix(NA, nrow=length(new.sp.names), ncol=length(new.ta.names),
		dimnames=list(new.sp.names, new.ta.names))

	# "merge". Overwrite with x information.
	newx[rownames(y), colnames(y)] <- y
	newx[rownames(x), colnames(x)] <- x

	# Create the new otuTable object
	newx <- taxTab(newx)

	return(newx)
})
################################################################################
#' @aliases merge_phyloseq_pair,sampleData,sampleData-method
#' @rdname merge_phyloseq_pair-methods
setMethod("merge_phyloseq_pair", signature("sampleData", "sampleData"), function(x, y){
	new.sa.names <- union(rownames(x), rownames(y))
	new.va.names <- union(colnames(x), colnames(y))
	
	partx <- data.frame("X0"=rownames(x), x)
	party <- data.frame("X0"=rownames(y), y)
	newx <- merge(partx, party, all=TRUE)
	# now we have the correct template, lets remove redundant rows.
	keep.samp.rows <- sapply(unique(as.character(newx[,1])), function(i,nx){
		rownames(subset(nx, X0==i))[1]
	},newx)
	newx <- newx[keep.samp.rows,]
	rownames(newx) <- as.character(newx$"X0")

	# "merge". Overwrite with x information.
	newx[rownames(y), colnames(y)] <- data.frame(y)
	newx[rownames(x), colnames(x)] <- data.frame(x)

	# trim the sample name column
	newx <- newx[,names(newx)!="X0"]
	
	# Create the new otuTable object
	newx <- sampleData(newx)
	return(newx)	
})
################################################################################
#' @aliases merge_phyloseq_pair,phylo,phylo-method
#' @rdname merge_phyloseq_pair-methods
setMethod("merge_phyloseq_pair", signature("phylo", "phylo"), function(x, y){
	ape::consensus(x, y)
})
################################################################################
################################################################################
#' Merge a subset of the species in \code{x} into one species/taxa/OTU.
#'
#' Takes as input an object that describes species/taxa
#' (e.g. \code{\link{phyloseq-class}}, \code{\link{otuTable-class}}, 
#'  \code{\link{phylo-class}}, \code{\link{taxonomyTable-class}}),
#' as well as 
#' a vector of species that should be merged.
#' It is intended to be able to operate at a low-level such that 
#' related methods, such as \code{\link{tipglom}} and \code{\link{taxglom}}
#' can both reliably call \code{merge_species} for their respective purposes.
#'
#' @usage merge_species(x, eqspecies, archetype=1)
#'
#' @param x (Required). An object that describes species (taxa). This includes
#'  \code{\link{phyloseq-class}}, \code{\link{otuTable-class}}, \code{\link{taxonomyTable-class}}, 
#'  \code{\link[ape]{phylo}}.
#' 
#' @param eqspecies (Required). The species names, or indices, that should be merged together.
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
#' @seealso \code{\link{tipglom}}, \code{\link{taxglom}}, \code{\link{merge_phyloseq}},
#'  \code{\link{merge_samples}}
#'
#' @export
#' @docType methods
#' @rdname merge_species-methods
#' @examples #
#' # # data(phylocom)
#' # # tree <- phylocom$phylo
#' # # otu  <- otuTable(phylocom$sample, speciesAreRows=FALSE)
#' # # otutree0 <- phyloseq(otu, tree)
#' # # plot(otutree0)
#' # # otutree1 <- merge_species(otutree0, tree$tip.label[1:8], 2)
#' # # plot(otutree1)
setGeneric("merge_species", function(x, eqspecies, archetype=1) standardGeneric("merge_species"))
###############################################################################
#' @aliases merge_species,otuTable-method
#' @rdname merge_species-methods
setMethod("merge_species", "otuTable", function(x, eqspecies, archetype=1){
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
# require(ape)
#' @aliases merge_species,phylo-method
#' @rdname merge_species-methods
setMethod("merge_species", "phylo", function(x, eqspecies, archetype=1){
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
################################################################################
#' @aliases merge_species,phyloseq-method
#' @rdname merge_species-methods
setMethod("merge_species", "phyloseq", function(x, eqspecies, archetype=1){
	comp_list   <- splat.phyloseq.objects(x)
	merged_list <- lapply(comp_list, merge_species, eqspecies, archetype)
	# the element names can wreak havoc on do.call
	names(merged_list) <- NULL
	# Re-instantiate the combined object using the species-merged object.
	do.call("phyloseq", merged_list)
})
###############################################################################
#' @aliases merge_species,sampleData-method
#' @rdname merge_species-methods
setMethod("merge_species", "sampleData", function(x, eqspecies, archetype=1){
	return(x)
})
###############################################################################
#' @aliases merge_species,taxonomyTable-method
#' @rdname merge_species-methods
setMethod("merge_species", "taxonomyTable", function(x, eqspecies, archetype=1){
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
# merge_species(ex4, species.names(ex4)[1:5])
# merge_species(phyloseqTree(ex4), species.names(ex4)[1:5])
# merge_species(phyloseq(ex4), 1:5)
# ################################################################################
# # Example of otuTree species merge
# otutree  = otuTree(ex4)
# otutree1 = merge_species(otutree, tre(otutree)$tip.label[1:9300])
# plot(tre(otutree1))
# # Not run, species indices not equivalent between phylo and otuTable.
# # Must use names (character):
# otutree2 = merge_species(otutree, 1:9300)
################################################################################
# # Examples of tree species merge:
# tree = tre(ex4)
# tree1 = merge_species(tree, tree$tip.label[1:500], 2)
# tree2 = merge_species(tree, 12:15, 1)
# tree3 = merge_species(tree, 12:15, 2)
# # Not run, won't know what merged:
# #tree3 = merge_species(tree, sample(1:32,15), 2)
# par(mfcol=c(2,2))
# plot(tree,  main="Tree 0")
# plot(tree1, main="Tree 1")
# plot(tree2, main="Tree 2")
# plot(tree3, main="Tree 3")
################################################################################
# # Example of otuTable species merge
# x4 = otuTable(matrix(sample(0:15,100,TRUE),40,10), speciesAreRows=TRUE)
# merge_species(x4, c("sp1", "sp3", "sp10", "sp28"), "sp10")
# merge_species(x4, c("sp1", "sp3", "sp10", "sp28", "sp35") )
# merge_species(x4, c("sp1", "sp3", "sp10", "sp28", "sp35") )@.Data
# merge_species(x4, 5:25)@.Data
# Not run:
# merge_species(x4, c("sp1", "sp3", "sp10", "sp28", "sp35", "") )
################################################################################
################################################################################
#' Merge samples based on a sample variable or factor.
#'
#' The purpose of this method is to merge/agglomerate the sample indices of a 
#' phyloseq object according to a categorical variable contained in a sampleData
#' or a provided factor.
#' 
#' NOTE: (\code{\link[ape]{phylo}}) trees and \code{\link{taxonomyTable-class}}
#' are not modified by this function, but returned in the output object as-is. 
#'
#' @usage merge_samples(x, group, fun=mean) 
#'
#' @param x (Required). An instance of a phyloseq class that has sample indices. This includes 
#'  \code{\link{sampleData-class}}, \code{\link{otuTable-class}}, and \code{\link{phyloseq-class}}. 
#'
#' @param group (Required). Either the a single character string matching a variable name in
#'  the corresponding sampleData of \code{x}, or a factor with the same length as
#'  the number of samples in \code{x}.
#'
#' @param fun (Optional). The function that will be used to merge the values that
#'  correspond to the same group for each variable. It must take a numeric vector
#'  as first argument and return a single value. Default is \code{\link[base]{mean}}.
#'  Note that this is (currently) ignored for the otuTable, where the equivalent
#'  function is \code{\link[base]{sum}}, but evaluated via \code{\link[base]{rowsum}}
#'  for efficiency.
#'
#' @return A phyloseq object that has had its sample indices merged according to
#'  the factor indicated by the \code{group} argument. The output class
#'  matches \code{x}.  
#'
#' @seealso \code{\link{merge_species}}, code{\link{merge_phyloseq}}
#'
#' @rdname merge_samples-methods
#' @docType methods
#' @export
#'
#' @examples #
#' # # # data(ex1)
#' # # # t1 <- merge_samples(sampleData(ex1), "Gender")
#' # # # t4 <- merge_samples(ex1, "Gender")
#' # # # identical(t1, sampleData(t4))
#' # # # t0 <- merge_samples(ex1, "Diet")
#' # # # grp <- as(data.frame(sampleData(ex1))[, "Diet"], "vector")
#' # # # t4 <- merge_samples(ex1, grp)
#' # # # identical(t0, t4)
#' # # # t1 <- merge_samples(otuTable(ex1), grp)
#' # # # t2 <- merge_samples(otuTable(ex1), factor(grp))
#' # # # t3 <- merge_samples(ex1, "Diet")
#' # # # identical(t1, t2)
#' # # # identical(t1, t3)
#' # # # identical(t1, otuTable(t3))
setGeneric("merge_samples", function(x, group, fun=mean) standardGeneric("merge_samples"))
################################################################################
#' @aliases merge_samples,sampleData-method
#' @rdname merge_samples-methods
setMethod("merge_samples", signature("sampleData"), function(x, group, fun=mean){
	x1    <- data.frame(x)

	# Check class of group and modify if "character"
	if( class(group)=="character" & length(group)==1 ){
		if( !group %in% colnames(x) ){stop("group not found among sample variable names.")}
		group <- x1[, group]
	}
	if( class(group)!="factor" ){
		# attempt to coerce to factor
		group <- factor(group)
	}

	# Remove any non-coercable columns.
	# Philosophy is to keep as much as possible. If it is coercable at all, keep.
	# Coerce all columns to numeric matrix
	coercable    <- sapply(x1, canCoerce, "numeric")
	x2           <- sapply(x1[, coercable], as, "numeric")
	rownames(x2) <- rownames(x1)	
	
	# Perform the aggregation.
	outdf <- aggregate(x2, list(group), fun)
	# get rownames from the "group" column (always first)
	# rownames(outdf) <- as.character(outdf[, 1])
	rownames(outdf) <- levels(group)
	# "pop" the first column
	outdf <- outdf[, -1, drop=FALSE]

	return( sampleData(outdf) )
})
################################################################################
#' @aliases merge_samples,otuTable-method
#' @rdname merge_samples-methods
setMethod("merge_samples", signature("otuTable"), function(x, group){
	# needs to be in sample-by-species orientation
	if( speciesAreRows(x) ){ x <- t(x) }
	# coerce to matrix, x2
	x2 <- as(x, "matrix")
	
	# # # #aggregate(x2, list(group), fun)
	out <- rowsum(x2, group)
	
	# convert back to otuTable, and return
	return( otuTable(out, speciesAreRows=FALSE) )
})
################################################################################
#' @aliases merge_samples,phyloseq-method
#' @rdname merge_samples-methods
setMethod("merge_samples", signature("phyloseq"), function(x, group, fun=mean){

	# Check if phyloseq object has a sampleData
	if( !is.null(access(x, "sampleData")) ){
		# Check class of group and modify if single "character" (column name)
		if( class(group)=="character" & length(group)==1 ){
			x1 <- data.frame(sampleData(x))		
			if( !group %in% colnames(x1) ){stop("group not found among sample variable names.")}
			group <- x1[, group]
		}
		if( class(group)!="factor" ){
			# attempt to coerce to factor
			group <- factor(group)
		}
		newSM <- merge_samples(sampleData(x), group, fun)
		newOT <- merge_samples(otuTable(x), group)
		phyloseqList <- list(newOT, newSM)
		
	# Else, the only relevant object to "merge_samples" is the otuTable
	} else {
		if( class(group)!="factor" ){ group <- factor(group) }
		phyloseqList <- list( newOT=merge_samples(otuTable(x), group) )
	}
	
	### Add to build-call-list the remaining components, if present in x.
	### NULL is returned by accessor if object lacks requested component/slot.
	### Order of objects in list doesn't matter for phyloseq.
	### The list should not be named.
	if( !is.null(access(x, "taxTab")) ){ phyloseqList <- c(phyloseqList, list(taxTab(x))) }
	if( !is.null(access(x, "tre"))    ){ phyloseqList <- c(phyloseqList, list(tre(x))) }
	
	return( do.call("phyloseq", phyloseqList) )
})
################################################################################