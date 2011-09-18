################################################################################
#' Agglomerate nearby tips in tree.
#' 
#' tipglom determines those groups of tips in the tree that are closer than the
#'   threshold provided as \code{speciationMinLength}. That is, all
#' tips of the tree separated by a cophenetic distance smaller than 
#' \code{speciationMinLength} will be agglomerated using \code{mergespecies}.
#' Can be used to create non-trivial OTU table, by agglomerating nearby tips.
#'
#' @param tree An object of class \code{otuTree} or its superclasses, in which case
#' the OTU argument can be ommitted. If, alternatively, \code{tree} is a \code{phylo}
#' or \code{phylo4} object, then \code{OTU} is required. 
#'
#' @param OTU An otuTable object. Optional. Ignored if \code{tree} is an 
#'  \code{otuTree} object or 
#'  one of its superclasses. If \code{tree} is a \code{phylo} or \code{phylo4}
#'  object and \code{OTU} is provided, then return will be an \code{otuTree}
#'  object. 
#'
#' @param speciationMinLength The minimum branch length that separates taxa. All
#' tips of the tree separated by a cophenetic distance smaller than 
#' \code{speciationMinLength} will be agglomerated. Default is 0.02
#'
#' @return An object of class \code{otuTree} or its superclasses. Output class matches
#' the class of \code{tree}, unless it is a \code{phylo}/\code{phylo4} object, in
#' which case \code{tipglom} returns an \code{otuTree} object.
#'
#' @rdname tipglom-methods
#' @docType methods
#' @export
#'
#' @examples #
#' # # # data(phylocom)
#' # # # otu  <- otuTable(phylocom$sample, speciesAreRows=FALSE)
#' # # # tree <- as(phylocom$phylo, "phylo4")
#' # # # x1   <- phyloseq(otu, tree)
#' # # # print(x1)
#' # # # library("phylobase")
#' # # # plot(tre(x1))
#' # # # x2 <- tipglom(x1, speciationMinLength=2.1)
#' # # # plot(tre(x2))
setGeneric("tipglom", function(tree, OTU, speciationMinLength=0.02) standardGeneric("tipglom"))
#' @rdname tipglom-methods
#' @aliases tipglom,phylo,otuTable-method
setMethod("tipglom", signature("phylo", "otuTable"), function(tree, OTU, speciationMinLength=0.02){
	# dispatch as the combined-object (otuTree-class), auto coherence.
	tipglom( phyloseq(OTU, tree), speciationMinLength=speciationMinLength)
})
#' @rdname tipglom-methods
#' @aliases tipglom,otuTree,ANY-method
setMethod("tipglom", signature("otuTree"), function(tree, speciationMinLength=0.02){
	tipglom.internal(tree, speciationMinLength=speciationMinLength)
})
#' @rdname tipglom-methods
#' @aliases tipglom,phylo,ANY-method
setMethod("tipglom", signature("phylo"), function(tree, speciationMinLength=0.02){
	tipglom.internal(tree, speciationMinLength=speciationMinLength)
})
#' @rdname tipglom-methods
#' @aliases tipglom,phylo4,ANY-method
setMethod("tipglom", signature("phylo4"), function(tree, speciationMinLength=0.02){
	tipglom.internal(tree, speciationMinLength=speciationMinLength)
})
################################################################################
#' Internal function for tiplgom.
#' 
#' Internal function, users should use the S4 method \code{\link{tipglom}}.
#' Tree can be higher-order object that contains a phylogenetic tree, 
#'   e.g. otuTree, etc. This is because \code{\link{mergespecies}} can
#' handle all the relevant objects, as can \code{\link{getTipDistMatrix}}.
#' Create Non-trivial OTU table, by agglomerating nearby tips.
#' tipglom.internal is called by the S4 \code{tipglom} methods. It is useful if 
#' a motivated user wants to see the internals of the implementation. By design
#' it lacks explicit object handling. Use \code{\link{tipglom}} instead.
#'
#' @param tree An object of class \code{phylo}, \code{phylo4}, \code{otuTree} 
#' or its superclasses. 
#'
#' @param speciationMinLength The minimum branch length that separates taxa. All
#' tips of the tree separated by a cophenetic distance smaller than 
#' \code{speciationMinLength} will be agglomerated.
#'
#' @return An object of class \code{phylo}, \code{phylo4}, \code{otuTree} 
#' or its superclasses. Output class matches the class of \code{tree}.
#'
#' @seealso tipglom
#' @import igraph
#'
#' @examples #
tipglom.internal = function(tree, speciationMinLength){
	# Create adjacency matrix, where tips are adjacent
	# if their distance is below the threshold speciationMinLength
	tipAdjacent <- (getTipDistMatrix( tree ) < speciationMinLength)
	# Define igraph object based on the tree-tip adjacenecy matrix
	ig          <- graph.adjacency(tipAdjacent, diag=FALSE)
	# Define the species cliques to loop through
	spCliques   <- edgelist2clique( get.edgelist(ig) )
	# successively merge
	for( i in 1:length(spCliques)){
		tree <- mergespecies(tree, eqspecies=spCliques[[i]])
	}
	# Test if you missed anything:
	# graph.adjacency( getTipDistMatrix(tree) < speciationMinLength, diag=FALSE )
	return(tree)
}
#################################################################
#' An internal wrapper function on \code{ape::cophenetic.phylo}
#' 
#' This is useful for determining tips that should be combined.
#' 
#' @param tree \code{phylo}
#' 
#' @param byRootFraction Should the distance be calculated according to
#' fractional distance to the root? If \code{FALSE}, the distance is
#' instead the patristic distance as calculated by cophenetic.phylo. 
#' Default \code{FALSE}.
#' 
#' @return character matrix. First column is the complete match, followed by
#'   one for each capture group
#' 
#' @seealso tipglom
#' @keywords internal
#' @aliases gettipdistmatrix getTipDistMatrix
#' @examples #
setGeneric("getTipDistMatrix", function(tree, byRootFraction=FALSE) standardGeneric("getTipDistMatrix"))
setMethod("getTipDistMatrix", signature("phylo"), function(tree, byRootFraction=FALSE){
	# require("picante")
	pairwiseSpecDists = cophenetic(tree)
	# If byRootFraction is true, normalize the cophenetic distances
	# according to the mean root age.
	if( byRootFraction ){
		# Want to normalize pairwise tip distances by their mean distance to root
		# start with tipAges
		tipAges = node.age(tree)$ages[which(tree$edge[,2] %in% 1:length(tree$tip.label))]
		names(tipAges) = tree$tip.label
		###### Want Mmean to be a matrix of the mean pairwise root-distance b/w each tip-pair
		Mmean = matrix(NA,length(tipAges),length(tipAges),
			dimnames=list(names(tipAges),names(tipAges))) 	
		means = combn(tipAges,2,mean)
		ind = combn(length(tipAges),2)
		for(i in 1:ncol(ind)){Mmean[ind[1,i], ind[2,i]] <- means[i]}
		for(i in 1:ncol(ind)){Mmean[ind[2,i], ind[1,i]] <- means[i]}
		diag(Mmean) <- tipAges
		# take the ratio of spec distances to the mean
		fracDists = pairwiseSpecDists / Mmean
		return(fracDists)
	} else {
		return(pairwiseSpecDists)
	}
})
setMethod("getTipDistMatrix", signature("phylo4"),
    function(tree, byRootFraction=FALSE){
    # Hand off to the "phylo" dispatch, rather than re-write for "phylo4"
	getTipDistMatrix(as(tree, "phylo"))
})
setMethod("getTipDistMatrix", signature("otuTree"), function(tree, byRootFraction=FALSE){
	getTipDistMatrix( tre(tree) )
})
gettipdistmatrix <- getTipDistMatrix
################################################################################
