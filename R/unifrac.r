##############################################################################
#' Internal function for \code{\link{wUniFrac}}.
#' 
#' A modified version of the \code{\link{internal2tips}} function, 
#' such that when a
#' tip is provided as int.node, that tip is returned. This is a more intuitive
#' behavior than the original picante version, which returns NULL.
#' This is currently used in \code{\link{wUniFrac}}.
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
	#require(picante); require(ape)	
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
    ### pjm2 add: If int.node was already a tip, then return itself, not NULL 
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
#' @keywords internal
#' @return character vector
#' @seealso wUniFrac
#' 
#' @examples #
UFwi = function(edge, samples, OTU, tree, AT=sum(OTU[samples[1],]), BT=sum(OTU[samples[2],])){
	A <- samples[1]
	B <- samples[2]
	# get the associated node in order to call internal2tips(). node number, n:
	n <- tree$edge[edge, 2]
	# get subset of tips (species) that descend from this branch
	edgeDescendants <- phyloseq:::internal2tips.self(tree, n, return.names = TRUE)
	# Calculate Ai/AT and Bi/BT, the number of individuals in each sample descended from branch n
	AiAT <- sum(OTU[A, edgeDescendants]) / AT
	BiBT <- sum(OTU[B, edgeDescendants]) / BT
	wi   <- abs(AiAT - BiBT)
	return(wi)	
}
##############################################################################
#' Calculate weighted UniFrac for a pair of samples from an OTU table.
#' 
#' Somewhat an internal function for \code{\link{wUniFrac}}, as usually interested
#' in the weighted-UniFrac distances between many different samples.
#' Returns a single numeric value, between 0 and 1 if normalized
#'
#' @usage wUniFracPair(OTU, tree, A, B, normalized=TRUE)
#'
#' @param OTU otuTable in samples-by-species orientation
#' @param tree object of class \code{phylo}
#' @param A single character string matching the first sample ID in the pair
#' @param B single character string matching the second sample ID in the pair
#' @param normalized Logical. Should the output be normalized such that values 
#'  range from 0 to 1 independent of branch length values? Default is \code{TRUE}. 
#'
#' @return numeric
#' @seealso wUniFrac
#' 
#' @export
#' @examples #
wUniFracPair = function(OTU, tree, A, B, normalized=TRUE){
	#require(picante); require(ape)
	# outer loop
	# Get the OTUs per sample. This should be calculated early so not needlessly repeated.
	AT <- sum(OTU[A,])
	BT <- sum(OTU[B,])
	# want to loop over each edge:
	edges <- 1:nrow(tree$edge)
	# sum the product of the branch length and their weights to get the numerator
	# numerator, u = sum( bi * UFwi ) = sum( bi * abs( Ai/AT - Bi/BT ) )
	numerator <- sum(tree$edge.length * sapply(edges, phyloseq:::UFwi, c(A,B), OTU, tree, AT, BT))
	# if not-normalized weighted UniFrac, just return "numerator";
	# the u-value in the w-UniFrac description
	if(!normalized){
		return(numerator)
	}
	if(normalized){
		# For denominator, we need the age of each tip
		# Get the tip ages from their associated edges (node.age gives the age of edges, ironically)
		tipAges <- picante::node.age(tree)$ages[which(tree$edge[, 2] %in% 1:length(tree$tip.label))]
		names(tipAges) <- tree$tip.label
		# denominator
		denominator <- sum( tipAges * (OTU[A, ]/AT + OTU[B, ]/BT) )
		# return the normalized weighted UniFrac values
		return(numerator / denominator)
	}
}
##############################################################################
##############################################################################
#' Calculate weighted UniFrac for all samples in an OTU table.
#'
#' This function calculates the weighted-UniFrac distance for all sample-pairs
#' in a species-abundance table using the abundances and a phylogenetic tree.
#' If \code{OTU} is a more complex object that already contains a phylogenetic
#' tree and abundance table, then the argument \code{tree} is not necessary
#' and will be ignored.  
#'
#' @usage wUniFrac(OTU, tree, normalized=TRUE, parallel=FALSE)
#'
#' @param OTU (Required). \code{otuTable}, or an \code{otuTree}. If you have
#'  instead a simple matrix of abundances, see \code{\link{otuTable}} for coercing
#'  it to the \code{otuTable} class.
#'
#' @param tree (Optional). Object of class \code{phylo}. Not necessary (and ignored) 
#'  if 
#'  \code{OTU} is an object that also contains a tree (\code{otuTree} class, 
#'  or its children).
#'
#' @param normalized (Optional). Logical. Should the output be normalized such that values 
#'  range from 0 to 1 independent of branch length values? Default is \code{TRUE}.
#'
#' @param parallel (Optional). Logical. Should execute calculation in parallel,
#'  using multiple CPU cores simultaneously? This can dramatically hasten the
#'  computation time for this function. However, it also requires that the user
#'  has registered a parallel ``backend'' prior to calling this function. 
#'  Default is \code{FALSE}. If FALSE, wUniFrac will register a serial backend
#'  so that \code{foreach::\%dopar\%} does not throw a warning.
#' 
#' @return a sample-by-sample distance matrix, suitable for NMDS, etc.
#' @seealso UniFrac unifrac vegdist
#' 
#' @docType methods
#' @export
#' @import foreach
#' @rdname wUniFrac-methods
#' @examples #
setGeneric("wUniFrac", function(OTU, tree, normalized=TRUE, parallel=FALSE) standardGeneric("wUniFrac"))
################################################################################
#' @aliases wUniFrac,otuTable,phylo-method
#' @rdname wUniFrac-methods
setMethod("wUniFrac", signature("otuTable", "phylo"),
										function(OTU, tree, normalized=TRUE, parallel=FALSE){
	#require(picante); require(ape)
    if (is.null(tree$edge.length)) {
        stop("Tree has no branch lengths, cannot compute UniFrac")
    }
    if (!is.rooted(tree)) {
        stop("Rooted phylogeny required for UniFrac calculation")
    }

	### Some parallel-foreach housekeeping.    
    # If user specifies not-parallel run (the default), register the sequential "back-end"
    if( !parallel ){ registerDoSEQ() }
    
	# create N x 2 matrix of all pairwise combinations of samples.
    spn <- combn(sample.names(OTU), 2, simplify=FALSE)
    # initialize UniFracMat with NAs
    UniFracMat <- matrix(NA, nsamples(OTU), nsamples(OTU))
    # define the rows/cols of UniFracMat with the sample names (rownames)    
    rownames(UniFracMat) <- sample.names(OTU)
    colnames(UniFracMat) <- sample.names(OTU)
    
	## Format coercion
	# Coerce to the picante orientation, with species as columns
	if( speciesAreRows(OTU) ){ OTU <- t( otuTable(OTU) ) }
   	# Coerce OTU to matrix for calculations.
    OTU <- as(OTU, "matrix")
	
   	# optionally-parallel implementation with foreach
  	distlist <- foreach( i = spn, .packages="phyloseq") %dopar% {
		A <- i[1]
		B <- i[2]
		return( phyloseq::wUniFracPair(OTU, tree, A, B, normalized) )
	}
	
    # This is in serial, but it is quick.
    distlist2distmat <- function(i, spn, DL){
    	UniFracMat[ spn[[i]][2], spn[[i]][1] ] <<- DL[[i]]
    }
	junk <- sapply(1:length(spn), distlist2distmat, spn, distlist)
    
    return(as.dist(UniFracMat))
})
################################################################################
#' @aliases wUniFrac,otuTree,ANY-method
#' @rdname wUniFrac-methods
#' @import phylobase
setMethod("wUniFrac", signature("otuTree"), 
		function(OTU, tree=NULL, normalized=TRUE, parallel=FALSE){
	# splat and pass to core function.
	wUniFrac(OTU=otuTable(OTU), tree=suppressWarnings(as(tre(OTU), "phylo")), normalized, parallel)
})
##############################################################################
##############################################################################
#' Calculate unweighted UniFrac from an abundance matrix and a tree.
#'
#' These methods depend on the picante and ape packages, as currently implemented.
#'
#' @usage UniFrac(OTU, tree)
#'
#' @param OTU otuTable, or a more-complex object that contains an otuTable.
#'
#' @param tree object of class \code{phylo}, defined by ape package.
#'
#' @return a distance matrix
#' @seealso wUniFrac unifrac
#' 
#' @docType methods
#' @rdname UniFrac-methods
#' @export
#' @examples #
setGeneric("UniFrac", function(OTU, tree) standardGeneric("UniFrac"))
###############################################################################
#' @aliases UniFrac,otuTable,phylo-method
#' @rdname UniFrac-methods
setMethod("UniFrac", signature("otuTable", "phylo"), function(OTU, tree){
	#require(picante)
    if (is.null(tree$edge.length)) {
        stop("Tree has no branch lengths, cannot compute UniFrac")
    }
    if (!is.rooted(tree)) {
        stop("Rooted phylogeny required for UniFrac calculation")
    }
    
    # Force orientation to "speciesAreCols"
    if( speciesarerows(OTU) ){OTU <- t(OTU)}
    
    # ensure that OTU is a matrix class OTU table
    OTU <- as(OTU, "matrix")
    
    # s is the number of samples. 
    # In Picante, the OTU-table is samples by species instead of species x samples
    s   <- nrow(OTU)
	# create N x 2 matrix of all pairwise combinations of samples.
    spn <- combn(rownames(OTU), 2)
    # OTU_comb sums the species abundances of all pair-wise combinations of samples from OTU
	# This is a way of combining presence/absence. Effictively "OR" logic for unweighted UniFrac
	OTU_comb <- t(apply(spn,2,function(i,OTU){ OTU[i[1],] + OTU[i[2],] },OTU))
    rownames(OTU_comb) <- apply(spn, 2, paste, sep="", collapse="_") 	
	# pdOTU is the total branch lengths of each sample
	pdOTU <- picante::pd(OTU, tree)
    # define pdOTU_comb, the total branch lengths for all combined sample pairs.
    pdOTU_comb <- picante::pd(OTU_comb, tree)
    # initialize UniFrac with NAs, a symmetric distance matrix
    UniFracMat <- matrix(NA, s, s)
    # define the rows/cols of UniFrac with the sample names (rownames)    
    rownames(UniFracMat) <- rownames(OTU)
    colnames(UniFracMat) <- rownames(OTU)  
	# Calculate unifrac at each index.
    for(i in 1:ncol(spn)){
		UniFracMat[spn[2,i], spn[1,i]] <- 2 - ( pdOTU[spn[1, i], "PD"] + 
			pdOTU[spn[2,i],"PD"] ) / pdOTU_comb[i, "PD"]
    }
    return(as.dist(UniFracMat))
})
###############################################################################
#' @aliases UniFrac,otuTree,ANY-method
#' @rdname UniFrac-methods
setMethod("UniFrac", signature("otuTree"), function(OTU, tree=NULL){
	UniFrac(otuTable(OTU), suppressWarnings(as(tre(OTU), "phylo")) )
})
################################################################################
