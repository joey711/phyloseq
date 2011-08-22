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
#' @param OTU optional. Ignored if \code{tree} is an \code{otuTree} object or 
#' one of its superclasses. Required if \code{tree} is a \code{phylo} or \code{phylo4}
#'  object. \code{OTU} must be an \code{otuTable} object or a 2-column table that
#' indicates the sample from which each sequence originated.
#'
#' @param speciationMinLength The minimum branch length that separates taxa. All
#' tips of the tree separated by a cophenetic distance smaller than 
#' \code{speciationMinLength} will be agglomerated. Default is 0.02
#'
#' @return An object of class \code{otuTree} or its superclasses. Output class matches
#' the class of \code{tree}, unless it is a \code{phylo}/\code{phylo4} object, in
#' which case \code{tipglom} returns an \code{otuTree} object.
#'
#' @export
#' @examples #
setGeneric("tipglom", function(tree, OTU, speciationMinLength=0.02) standardGeneric("tipglom"))
setMethod("tipglom", signature("phylo", "otuTable"), function(tree, OTU, speciationMinLength=0.02){
	# pass off to the combined-object (otuTree) version
	tipglom( phyloseq(OTU, tree), speciationMinLength=speciationMinLength)
})
setMethod("tipglom", signature("phylo4", "otuTable"), function(tree, OTU, speciationMinLength=0.02){
  # pass off to the combined-object (otuTree) version
	tipglom( phyloseq(OTU, tree), speciationMinLength=speciationMinLength)
})
setMethod("tipglom", signature("otuTree"), function(tree, speciationMinLength=0.02){
	tipglom.internal(tree, speciationMinLength=speciationMinLength)
}) 
setMethod("tipglom", signature("phylo"), function(tree, speciationMinLength=0.02){
	tipglom.internal(tree, speciationMinLength=speciationMinLength)
})
setMethod("tipglom", signature("otuTree4"), function(tree, speciationMinLength=0.02){
	tipglom.internal(tree, speciationMinLength=speciationMinLength)
})
setMethod("tipglom", signature("phylo4"), function(tree, speciationMinLength=0.02){
	tipglom.internal(tree, speciationMinLength=speciationMinLength)
})
################################################################################
#' Internal function for tiplgom.
#' 
#' Internal function, users should use the S4 method \code{\link{tipglom}}.
#' Tree can be higher-order object that contains a phylogenetic tree, 
#'   e.g. otuTree, otuTree4, etc. This is because \code{\link{mergespecies}} can
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
#' A wrapper function on cophenetic.phylo()
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
setMethod("getTipDistMatrix", signature("otuTree4"), function(tree, byRootFraction=FALSE){
  getTipDistMatrix( tre(tree) )
})
gettipdistmatrix <- getTipDistMatrix
################################################################################
# Examples:
# plot(ig, layout=layout.fruchterman.reingold, 
	# vertex.size=0.6, vertex.label.dist=0.1, 
	# edge.arrow.mode="-",vertex.color="red",
	# vertex.label=NA,edge.color="blue")
################################################################################
# Test igraph2cliques later. Not sure if that hard indexing works on igraph
# objects. Values don't match edgelist2clique
# ################################################################################
# igraph2cliques <- function(ig){
	# # ig is a igraph "graph" object
	# require("statnet")
	# igfac = factor(label.propagation.community(ig))
	# spCliques = ig[[9]][[3]][[1]][igfac]
	# names(spCliques) = ig[[9]][[3]][[1]]
	# tabCliques = table(spCliques) 
	# tabCliques = tabCliques[tabCliques > 1]
	# cliqueNames = names(factor(tabCliques))
	# cliques = lapply(cliqueNames, function(i,spCliq){
		# names(spCliques[spCliq %in% i])
	# }, spCliques)
	# return(cliques)
# }
# cl1 = igraph2cliques(ig)
################################################################################
