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
#' ## map1 <- sampleMap(map1)
#' ## ex1 <- phyloseq(OTU1, map1, tax1)
#' ## x <- ex1
#' ## x <- phyloseq(ex1)
#' ## y <- taxTab(ex1)
#' ## merge_phyloseq(x, y)
#' ## merge_phyloseq(y, y, y, y)
#' ## # Try the simple example dataset, ex1
#' ## data(ex1)
#' ## merge_phyloseq( otutree(ex1), sampleMap(ex1))
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
	# Remove names to avoid any conflicts with phyloseq(), which does not have named-arguments
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
#' Special note: trees are merged using \code{ape::\link{consensus}}.
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
#' @seealso \code{\link{merge_phyloseq}} \code{\link{mergespecies}}
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
#' ## # merge two simulated sampleMap objects
#' ## x <- data.frame( matrix(sample(0:3,250,TRUE),25,10), 
#' ##   matrix(sample(c("a","b","c"),150,TRUE),25,6) )
#' ## x <- sampleMap(x)
#' ## y <- data.frame( matrix(sample(4:6,200,TRUE),20,10), 
#' ##   matrix(sample(c("d","e","f"),120,TRUE),20,8) )
#' ## y <- sampleMap(y)
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
#' @aliases merge_phyloseq_pair,sampleMap,sampleMap-method
#' @rdname merge_phyloseq_pair-methods
setMethod("merge_phyloseq_pair", signature("sampleMap", "sampleMap"), function(x, y){
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
	newx <- sampleMap(newx)
	return(newx)	
})
################################################################################
#' @aliases merge_phyloseq_pair,phylo,phylo-method
#' @rdname merge_phyloseq_pair-methods
setMethod("merge_phyloseq_pair", signature("phylo", "phylo"), function(x, y){
	consensus(x, y)
})
################################################################################
#' @aliases merge_phyloseq_pair,phylo4,phylo4-method
#' @rdname merge_phyloseq_pair-methods
setMethod("merge_phyloseq_pair", signature("phylo4", "phylo4"), function(x, y){
	x <- as(x, "phylo")
	y <- as(y, "phylo")
	# dispatch to the method for "phylo" trees.
	phylo_tree <- merge_phyloseq_pair(x, y)
	# convert back to phylo4 and return.
	return( as(phylo_tree, "phylo4") )
})
################################################################################
