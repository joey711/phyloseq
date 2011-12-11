##############################################################################
#' Custom version of \code{picante::internal2tips()}
#'
#' Internal function for \code{\link{UniFrac}}.
#' 
#' A modified version of the \code{\link{internal2tips}} function, 
#' such that when a
#' tip is provided as \code{int.node}, that tip is returned. This is a more intuitive
#' behavior than the original picante version, which returns NULL.
#' This is currently used in \code{\link{UniFrac}}.
#'
#' @param phy object of class \code{phylo}.
#'
#' @param int.node the internal node you want to get the descendants of.
#'
#' @keywords internal
#' @return character vector
#' @seealso internal2tips UniFrac
#' 
#' @examples #
internal2tips.self = function (phy, int.node, return.names = TRUE){
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
#' Internal function for unweighted UniFrac edge-weight.
#' 
#' A function that takes a phylo object (tree) and an edge-index
#' as input, and returns the edge-weight term for unweighted UniFrac.
#'
#' @param edge The edge index
#'
#' @param samples Character vector of length 2, giving the pair of
#'  samples under comparison.
#'
#' @param occ \code{otuTable} object in samples-by-species orientation
#'
#' @param tree object of class \code{phylo}
#'
#' @keywords internal
#' @return character vector
#' @seealso UniFrac
#' 
#' @examples #
ufnum <- function(edge, samples, occ, tree){
	A <- samples[1]
	B <- samples[2]
	
	# get the associated node in order to call internal2tips(). node number, n:
	n <- tree$edge[edge, 2]
	
	# get subset of tips (species) that descend from this branch
	edgeDescendants <- phyloseq:::internal2tips.self(tree, n, return.names = TRUE)
	
	# Tabulate whether Ai and Bi have descendants in each sample descended from branch n
	Ai <- (sum(occ[A, edgeDescendants]) > 0) - 0
	Bi <- (sum(occ[B, edgeDescendants]) > 0) - 0
	wi   <- abs(Ai - Bi) > 0
	return(wi)	
}
################################################################################
#' Internal function for weighted UniFrac edge-weight.
#' 
#' A function that takes a phylo object (tree) and an edge-index
#' as input, and returns the edge-weight term for weighted \code{\link{UniFrac}}.
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
#' @seealso UniFrac
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
################################################################################
#' Calculate (unweighted) UniFrac for a pair of samples from an occurrence matrix.
#' 
#' Somewhat an internal function for \code{\link{UniFrac}}, as 
#' an investigator is usually interested
#' in the UniFrac distances between many different samples.
#' Returns a single numeric value, between 0 and 1.
#'
#' @usage unifracPair(occ, tree, A, B)
#'
#' @param occ (Required). A matrix. A samples-by-species occurrence matrix.
#' @param tree (Required). A rooted phylogenetic tree. \code{phylo} class.
#' @param A (Required). Single character string. The name of sample \code{"A"}.
#' @param B (Required). Single character string. The name of sample \code{"B"}. 
#'
#' @return A single number between 0, 1.
#' @seealso See the main function, \code{\link{UniFrac}}.
#' 
#' @export
#' @examples #
unifracPair <- function(occ, tree, A, B){
	
	# Prune tree to just those species present in either A or B
	ispecies <- colnames(occ)[ (occ[A, ] > 0) | (occ[B, ] > 0) ]
	itree    <- prune_species(ispecies, tree)
	
	# numerator loops over each edge:
	edges <- 1:nrow(itree$edge)
		
	# sum the abs value of the occurrence matrix
	# numerator, u = sum( bi * abs( Ai - Bi )  )
	#numerator <- sum(tree$edge.length * sapply(edges, phyloseq:::ufnum, c(A,B), occ, tree))
	numerator <- sum(itree$edge.length[ sapply(edges, phyloseq:::ufnum, c(A, B), occ, itree) ])
	
	# denominator is the sum of the branch lengths for which A and B have species
	# (not all branches. Branches having nothing to do with A,B must not count)
	return( numerator / sum(itree$edge.length) )
}
##############################################################################
#' Calculate weighted UniFrac for a pair of samples from an abundance matrix.
#' 
#' Somewhat an internal function for \code{\link{UniFrac}}, as usually interested
#' in the weighted-UniFrac distances between many different samples.
#' Returns a single numeric value, between 0 and 1 if normalized
#'
#' @usage wUniFracPair(OTU, tree, A, B, normalized=TRUE)
#'
#' @param OTU (Required). An abundance matrix in samples-by-species orientation.
#' @param tree (Required). Object of class \code{\link[ape]{phylo}}
#' @param A (Required). single character string matching the first sample ID in the pair
#' @param B (Required). single character string matching the second sample ID in the pair
#' @param normalized (Optional). Logical. Should the output be normalized such that values 
#'  range from 0 to 1 independent of branch length values? Default is \code{TRUE}. 
#'
#' @return A single number between 0, 1.
#' @seealso UniFrac
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
#' Calculate weighted or unweighted UniFrac for all samples in an OTU table.
#'
#' This function calculates the UniFrac distance for all sample-pairs
#' in a species-abundance table using the abundances and a phylogenetic tree.
#' If \code{OTU} is a more complex object that already contains a phylogenetic
#' tree and abundance table, then the argument \code{tree} is not necessary
#' and will be ignored. WARNING: The species names of both \code{OTU} and
#' \code{tree} must match exactly, or weird index erros can result. A warning
#' has been added to further protect users from encountering this issue unknowingly.
#' The easiest solution is to combine the \code{OTU} and \code{tree} objects 
#' using the \code{\link{phyloseq}} function. E.g. \code{phyloseq(OTU, tree)}
#' will return an \code{phyloseq}-class object that has been pruned and comprises
#' the minimum arguments necessary for \code{UniFrac()}. 
#'
#' Parallelization is possible, and encouraged for this computing-intensive calculation.
#' It has been implemented with the \emph{foreach} package.
#' This means that parallel calls need to be preceded by 2 or more commands that
#' that register the parallel ``backend''. This is acheived via your choice of
#' helper packages. One of the simplest seems to be the \emph{doMC} package.
#' See the commented-out examples at the bottom for running \code{UniFrac()}
#' in parallel via doMC. 
#'
#' For more information, see the following links on registering the ``backend'':
#' 
#' \emph{foreach} package manual: 
#'
#' \url{http://cran.r-project.org/web/packages/foreach/index.html}
#'
#' Notes on parallel computing in \code{R}. Skip to the section describing 
#' the \emph{foreach Framework}. It gives off-the-shelf examples for registering
#' a parallel backend using the \emph{doMC}, \emph{doSNOW}, or \emph{doMPI} packages:
#' 
#' \url{http://trg.apbionet.org/euasiagrid/docs/parallelR.notes.pdf}
#'
#' @usage UniFrac(OTU, tree, weighted=FALSE, normalized=TRUE, parallel=FALSE)
#'
#' @param OTU (Required). \code{otuTable}, or \code{phyloseq} object. If you have
#'  instead a simple matrix of abundances, see \code{\link{otuTable}} for coercing
#'  it to the \code{otuTable} class.
#'
#' @param tree (Optional). Object of class \code{phylo}. Not necessary (and ignored) 
#'  if 
#'  \code{OTU} is an object that also contains a tree (\code{phyloseq} class, 
#'  or its children).
#'
#' @param weighted (Optional). Logical. Should use weighted-UniFrac calculation?
#'  Weighted-UniFrac takes into account the relative abundance of species/taxa
#'  shared between samples, whereas unweighted-UniFrac only considers 
#'  presence/absence. Default is \code{FALSE}, meaning the unweighted-UniFrac
#'  distance is calculated for all pairs of samples.
#'
#' @param normalized (Optional). Logical. Should the output be normalized such that values 
#'  range from 0 to 1 independent of branch length values? Default is \code{TRUE}.
#'  Note that (unweighted) \code{UniFrac} is always normalized by total branch-length,
#'  and so this value is ignored when \code{weighted == FALSE}.
#'
#' @param parallel (Optional). Logical. Should execute calculation in parallel,
#'  using multiple CPU cores simultaneously? This can dramatically hasten the
#'  computation time for this function. However, it also requires that the user
#'  has registered a parallel ``backend'' prior to calling this function. 
#'  Default is \code{FALSE}. If FALSE, UniFrac will register a serial backend
#'  so that \code{foreach::\%dopar\%} does not throw a warning.
#' 
#' @return a sample-by-sample distance matrix, suitable for NMDS, etc.
#' @seealso \code{vegan::vegdist}
#'
#' @references \url{http://bmf.colorado.edu/unifrac/}
#'
#' Lozupone, Hamady and Knight, ``UniFrac - An Online Tool for Comparing Microbial 
#' Community Diversity in a Phylogenetic Context.'', BMC Bioinformatics 2006, 7:371
#'
#' Lozupone, Hamady, Kelley and Knight, ``Quantitative and qualitative (beta) 
#' diversity measures lead to different insights into factors that structure
#' microbial communities.'' Appl Environ Microbiol. 2007
#'
#' Lozupone C, Knight R. ``UniFrac: a new phylogenetic method for comparing microbial
#' communities.'' Appl Environ Microbiol. 2005 71 (12):8228-35.
#' 
#' @docType methods
#' @export
#' @import foreach
#' @rdname UniFrac-methods
#' @examples #
#' # # ########################################
#' # # # Create a simulated tree and abundance matrix
#' # # ########################################
#' # # ### Set dimensions of random data.
#' # # Nspec <- 50
#' # # Nsamp <- 10
#' # # ### Create random otuTable
#' # # abvec <- sample(c(rep(0, 5), 1:10), Nsamp*Nspec, TRUE)
#' # # rOTU  <- otuTable(matrix(abvec, nrow=Nsamp, ncol=Nspec), speciesAreRows=FALSE)
#' # # ### Create random tree
#' # # rtree <- rcoal(Nspec, tip.label=species.names(rOTU))
#' # # ### Combine into one object.
#' # # rotuTree <- phyloseq(rOTU, rtree)

#' # # ########################################
#' # # # Run UniFrac
#' # # ########################################
#' # # ### Unweighted UniFrac
#' # # unweightedUF.dist <- UniFrac(rotuTree)

#' # # ### weighted UniFrac
#' # # weightedUF.ser.dist <- UniFrac(rotuTree, weighted=TRUE)

#' # # ### parallel weighted UniFrac
#' # # require("doMC")
#' # # registerDoMC(4)
#' # # getDoParWorkers()
#' # # weightedUF.par.dist <- UniFrac(rotuTree, weighted=TRUE, parallel=TRUE)

#' # # ### Compare with picante version
#' # # picUF <- unifrac( as(rOTU, "matrix"), rtree)
#' # # round((unweightedUF.dist - picUF), 5)
#' ##
setGeneric("UniFrac", function(OTU, tree, weighted=FALSE, normalized=TRUE, parallel=FALSE){
	standardGeneric("UniFrac")
})
################################################################################
#' @aliases UniFrac,otuTable,phylo-method
#' @rdname UniFrac-methods
setMethod("UniFrac", signature("otuTable", "phylo"),
		function(OTU, tree, weighted=FALSE, normalized=TRUE, parallel=FALSE){
											
	# # # require(picante); require(ape); require(foreach)

	# Some important checks.
    if( is.null(tree$edge.length) ) {
        stop("Tree has no branch lengths, cannot compute UniFrac")
    }
    if( !is.rooted(tree) ) {
        stop("Rooted phylogeny required for UniFrac calculation")
    }
    # Attempt to prune if species-names do not match between tree and OTU.
    if( !setequal(species.names(tree), species.names(OTU)) ){
		warning("Species names between tree and OTU do not match exactly. They must.\n")
		warning("Attempting to prune non-shared species from either object... \n")
		warning("To solve, try: phyloseq(OTU, tree)    or try:  prune_species( ).\n")
		combinedOTUtree <- phyloseq(OTU, tree)
		OTU  <- otuTable(combinedOTUtree)
		tree <- suppressWarnings(as(tre(combinedOTUtree), "phylo"))
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

	if(!weighted){
		# Coerce OTU to occurrence matrix for unweighted-UF calculations
		occ <- (OTU > 0) - 0
	}
	
   	# optionally-parallel implementation with foreach
  	distlist <- foreach( i = spn, .packages="phyloseq") %dopar% {
		A <- i[1]
		B <- i[2]
		if( weighted ){
			return( phyloseq::wUniFracPair(OTU, tree, A, B, normalized) )
		} else {
			return( phyloseq::unifracPair(occ, tree, A, B) )
		}
	}
	
    # This is in serial, but it is quick.
    distlist2distmat <- function(i, spn, DL){
    	UniFracMat[ spn[[i]][2], spn[[i]][1] ] <<- DL[[i]]
    }
	junk <- sapply(1:length(spn), distlist2distmat, spn, distlist)
    
    return(as.dist(UniFracMat))
})
################################################################################
#' @aliases UniFrac,phyloseq,ANY-method
#' @rdname UniFrac-methods
setMethod("UniFrac", signature("phyloseq"), 
		function(OTU, tree=NULL, weighted=FALSE, normalized=TRUE, parallel=FALSE){
	# splat and pass to core function.
	UniFrac(OTU=otuTable(OTU), tree=tre(OTU), weighted, normalized, parallel)
})
##############################################################################
