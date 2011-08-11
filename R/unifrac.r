##############################################################################
#' Internal function for \code{\link{wUniFrac}}.
#' 
#' A modified version of the \code{\link{internal2tips}} function, 
#' such that when a
#' tip is provided as int.node, that tip is returned. This is a more intuitive
#' behavior than the picante version, which returns NULL, and currently used
#' in the \code{\link{wUniFrac}}.
#'
#' @param phy object of class \code{phylo}.
#'
#' @param int.node the internal node you want to get the descendants of.
#'
#' @keywords internal
#' @return character vector
#' @seealso internal2tips wUniFrac
#' 
#' @examples #
internal2tips.self = function (phy, int.node, return.names = FALSE){
	require(picante); require(ape)	
    Ntaxa = length(phy$tip.label)
    Nnode = phy$Nnode
    if ((Ntaxa + Nnode - 1) != nrow(phy$edge)) {
        print("tree structure error")
        break
    }
    # check if character index of node was used and convert, else use int.mode as index
    if (mode(int.node) == "character"){
		nodes = which(phy$node.label == int.node) + Ntaxa
	}else nodes = int.node
	
    # initialize the return character vector "tips"
    tips = c()
    
    # loop repeats while there are internal to include 
    repeat {
        nodes = phy$edge[which(phy$edge[, 1] %in% nodes), 2]
        if (length(nodes) == 0) 
            break
        tips = c(tips, nodes)
    }
    # Remove any internal nodes. They have larger indices than the number of taxa (tips)
    tips = tips[tips <= Ntaxa]
    # JOEY ADD: If int.node was already a tip, then return itself, not NULL 
    if( int.node <= Ntaxa & length(tips) == 0 ){
    	tips = int.node
    }
    if (return.names) 
        tips = phy$tip.label[tips]
    return(tips)
}
##############################################################################
#' Internal function for \code{wUniFrac}.
#' 
#' A function that takes a phylo object (tree) and an edge-index
#' as input, and returns the edge-weight term for weighted UniFrac.
#'
#' @param edge The edge index
#'
#' @param samples Character vector of length 2, giving the pair of
#'  samples under comparison.
#'
#' @param OTU \code{otuTable} object in samples-by-species orientation
#'
#' @param tree object of class \code{phylo}
#'
#' @keywords internal internal2tips
#' @return character vector
#' @seealso wUniFrac
#' 
#' @examples #
UFwi = function(edge,samples,OTU,tree,AT=sum(OTU[samples[1],]),BT=sum(OTU[samples[2],])){
	#edge = edges[1]
	A = samples[1]; B = samples[2]
	# get the associated node in order to call internal2tips(). node number, n:
	n = tree$edge[edge,2]
	# get subset of tips (species) that descend from this branch
	edgeDescendants = internal2tips.self(tree, n, return.names = TRUE)
	# Calculate Ai/AT and Bi/BT, the number of individuals in each sample descended from branch n
	AiAT = sum(OTU[A,edgeDescendants]) / AT
	BiBT = sum(OTU[B,edgeDescendants]) / BT
	wi = abs(AiAT - BiBT)
	return(wi)	
}
##############################################################################
#' Calculate weighted UniFrac for a pair of samples from an OTU table.
#' 
#' Somewhat an internal function for \code{\link{wUniFrac}}, as usually interested
#' in the weighted-UniFrac distances between many different samples.
#' Returns a single numeric value, between 0 and 1 if normalized
#'
#' @param OTU otuTable in samples-by-species orientation
#' @param tree object of class \code{phylo}
#' @param A single character string matching the first sample ID in the pair
#' @param B single character string matching the second sample ID in the pair
#' @param UFwi The UFwi function. May be unnecessary. Fixes namespace issues
#'  during certain apply loops. The function UFwi should be provided.
#' @param normalized Logical. Should the output be normalized such that values 
#'  range from 0 to 1 independent of branch length values? Default is \code{TRUE}. 
#'
#' @return numeric
#' @seealso wUniFrac
#' 
#' @export
#' @examples #
wUniFracPair = function(OTU, tree, A, B, UFwi=UFwi, normalized=TRUE){
	require(picante); require(ape)
	# outer loop
	# Get the OTUs per sample. This should be calculated early so not needlessly repeated.
	AT = sum(OTU[A,]); BT = sum(OTU[B,])
	# want to loop over each edge:
	edges = 1:nrow(tree$edge)
	# sum the product of the branch length and their weights to get the numerator
	# numerator, u = sum( bi * UFwi ) = sum( bi * abs( Ai/AT - Bi/BT ) )
	numerator = sum(tree$edge.length * sapply(edges,UFwi,c(A,B),OTU,tree,AT,BT))
	# if not-normalized weighted UniFrac, just return "numerator";
	# the u-value in the w-UniFrac description
	if(!normalized){
		return(numerator)
	}
	if(normalized){
		# For denominator, we need the age of each tip
		# Get the tip ages from their associated edges (node.age gives the age of edges, ironically)
		tipAges = node.age(tree)$ages[which(tree$edge[,2] %in% 1:length(tree$tip.label))]
		names(tipAges) = tree$tip.label
		# denominator
		denominator = sum( tipAges * (OTU[A,]/AT + OTU[B,]/BT) )
		# return the normalized weighted UniFrac values
		return(numerator / denominator)
	}
}
##############################################################################
#' Calculate weighted UniFrac for all samples in an OTU table.
#'
#' This function calculates the weighted-UniFrac distance for all sample-pairs
#' in a species-abundance table using the abundances and a phylogenetic tree.
#' If \code{OTU} is a more complex object that already contains a phylogenetic
#' tree and abundance table, then the argument \code{tree} is not necessary
#' and will be ignored.  
#'
#' @param OTU otuTable in samples-by-species orientation.
#'
#' @param tree object of class \code{phylo}
#'
#' @param normalized Logical. Should the output be normalized such that values 
#'  range from 0 to 1 independent of branch length values? Default is \code{TRUE}.
#' 
#' @keywords UniFrac weighted-unifrac
#' 
#' @return a sample-by-sample distance matrix, suitable for NMDS, etc.
#' @seealso UniFrac unifrac vegdist
#' 
#' @import picante ape
#' @export
#' @examples #
wUniFrac = function(OTU, tree, normalized=TRUE){
	#require(picante); require(ape)
    if (is.null(tree$edge.length)) {
        stop("Tree has no branch lengths, cannot compute UniFrac")
    }
    if (!is.rooted(tree)) {
        stop("Rooted phylogeny required for UniFrac calculation")
    }
	# ensure that OTU is a matrix class OTU table
    OTU <- as.matrix(OTU)
    # s is the number of samples. In Picante, the OTU-table is samples by species
    # instead of species by samples
    s <- nrow(OTU)
	# create N x 2 matrix of all pairwise combinations of samples.
    spn = combn(rownames(OTU),2)
    # initialize UniFrac with NAs, a symmetric distance matrix
    UniFrac <- matrix(NA, s, s)
    # define the rows/cols of UniFrac with the sample names (rownames)    
    rownames(UniFrac) <- rownames(OTU)
    colnames(UniFrac) <- rownames(OTU)
  	# Calculate unifrac at each index.
    for(i in 1:ncol(spn)){
		A = spn[1,i]; B = spn[2,i]
    	UniFrac[B,A] = wUniFracPair(OTU,tree,A,B,UFwi,normalized)
    }
    return(as.dist(UniFrac))
}
##############################################################################
#' Alternative 'vectorized' implementation of the (unweighted) UniFrac calculation
#'
#' @param OTU otuTable in samples-by-species orientation
#'
#' @param tree object of class \code{phylo}
#'
#' @keywords UniFrac
#' @return a distance matrix
#' @seealso wUniFrac unifrac
#' 
#' @export
#' @examples #
UniFrac = function(OTU, tree){
	require(picante)
    if (is.null(tree$edge.length)) {
        stop("Tree has no branch lengths, cannot compute UniFrac")
    }
    if (!is.rooted(tree)) {
        stop("Rooted phylogeny required for UniFrac calculation")
    }
    # ensure that OTU is a matrix class OTU table
    OTU <- as.matrix(OTU)
    # s is the number of samples. In Picante, the OTU-table is samples x species instead of species x samples
    s <- nrow(OTU)
	# create N x 2 matrix of all pairwise combinations of samples.
    spn = combn(rownames(OTU), 2)
    # OTU_comb sums the species abundances of all pair-wise combinations of samples from OTU
	# This is a way of combining presence/absence. Effictively "OR" logic for unweighted UniFrac
	OTU_comb = t(apply(spn,2,function(i,OTU){ OTU[i[1],] + OTU[i[2],] },OTU))
    rownames(OTU_comb) = apply(spn,2,paste,sep="",collapse="_") 	
	# pdOTU is the total branch lengths of each sample
	pdOTU <- pd(OTU, tree)
    # define pdOTU_comb, the total branch lengths for all combined sample pairs.
    pdOTU_comb <- pd(OTU_comb, tree)
    # initialize UniFrac with NAs, a symmetric distance matrix
    UniFrac <- matrix(NA, s, s)
    # define the rows/cols of UniFrac with the sample names (rownames)    
    rownames(UniFrac) <- rownames(OTU)
    colnames(UniFrac) <- rownames(OTU)  
	# Calculate unifrac at each index.
    for(i in 1:ncol(spn)){
    	UniFrac[spn[2,i], spn[1,i]] <- 2 - ( pdOTU[spn[1, i], "PD"] + 
        pdOTU[spn[2,i],"PD"] ) / pdOTU_comb[i, "PD"]
    }
    return(as.dist(UniFrac))
}
################################################################################
# Extend UniFrac methods to handle phyloseq classes. 
################################################################################
setGeneric("UniFrac", function(OTU, tree) standardGeneric("UniFrac"))
setMethod("UniFrac", signature("otuTable", "phylo"), function(OTU, tree){
  # not yet re-written. Needs to be parallelized.
  # Just borrow regular unifrac for now.
  require("picante")
  if ( speciesAreRows(OTU) ){
		OTU = t( OTU )
	}
	unifrac(OTU@.Data, tree)  
})
setMethod("UniFrac", signature("otuTree"), function(OTU, tree=NULL){
	UniFrac(otuTable(OTU), tre(OTU))
})
################################################################################
# Extend wUniFrac methods to handle phyloseq class. 
################################################################################
setGeneric("wUniFrac")
setMethod("wUniFrac", signature("otuTree"),
  function(OTU, tree=NULL, normalized=TRUE){
	if ( speciesAreRows(OTU) ){
		otu <- t( otuTable(OTU) )
	} else { otu <- otuTable(OTU) }
	wUniFrac(otu, tre(OTU), normalized)
})
################################################################################
