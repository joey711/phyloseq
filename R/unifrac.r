##############################################################################
#' Custom version of \code{\link[picante]{internal2tips}}
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
#' @seealso \code{\link[picante]{internal2tips}} \code{\link{UniFrac}}
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
ufnum <- function(edge, samples, occ, tree){
	A <- samples[1]
	B <- samples[2]
	
	# get the associated node in order to call internal2tips(). node number, n:
	n <- tree$edge[edge, 2]
	
	# get subset of tips (species) that descend from this branch
	edgeDescendants <- internal2tips.self(tree, n, return.names = TRUE)
	
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
UFwi = function(edge, samples, OTU, tree, AT=sum(OTU[samples[1],]), BT=sum(OTU[samples[2],])){
	A <- samples[1]
	B <- samples[2]
	# get the associated node in order to call internal2tips(). node number, n:
	n <- tree$edge[edge, 2]
	# get subset of tips (species) that descend from this branch
	edgeDescendants <- internal2tips.self(tree, n, return.names = TRUE)
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
#' @keywords internal
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
#' @seealso \code{\link{UniFrac}}
#' @keywords internal
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
		# denominator (assumes tree-indices and otuTable indices are same order)
		denominator <- sum( tipAges * (OTU[A, ]/AT + OTU[B, ]/BT) )
		# return the normalized weighted UniFrac values
		return(numerator / denominator)
	}
}
##############################################################################
#' @keywords internal
originalUniFrac <- function(physeq, weighted=FALSE, normalized=TRUE, parallel=FALSE, fast=TRUE){
	# # # require(picante); require(ape); require(foreach)
	# Access the needed components. Note, will error if missing in physeq.
	OTU  <- otuTable(physeq)
	tree <- tre(physeq)

	# Some important checks.
    if( is.null(tree$edge.length) ) {
        stop("Tree has no branch lengths, cannot compute UniFrac")
    }
    if( !is.rooted(tree) ) {
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
   	# Enforce that tree and otuTable indices are the same order, 
   	# by re-ordering OTU if needed
	if( !all(colnames(OTU) == species.names(tree)) ){
		OTU <- OTU[, species.names(tree)]
	}    
	# If unweighted-UniFrac, coerce to a presence-absence contingency, occ
	if(!weighted){
		# Coerce OTU to occurrence matrix for unweighted-UF calculations
		occ <- (OTU > 0) - 0
	}
	
   	# optionally-parallel implementation with foreach
  	distlist <- foreach( i = spn, .packages="phyloseq") %dopar% {
		A <- i[1]
		B <- i[2]
		if( weighted ){
			return( phyloseq:::wUniFracPair(OTU, tree, A, B, normalized) )
		} else {
			return( phyloseq:::unifracPair(occ, tree, A, B) )
		}
	}
	
    # This is in serial, but it is quick.
    distlist2distmat <- function(i, spn, DL){
    	UniFracMat[ spn[[i]][2], spn[[i]][1] ] <<- DL[[i]]
    }
	junk <- sapply(1:length(spn), distlist2distmat, spn, distlist)
    
    return(as.dist(UniFracMat))

}
##############################################################################
#' Calculate weighted or unweighted (Fast) UniFrac distance for all sample pairs.
#'
#' This function calculates the (Fast) UniFrac distance for all sample-pairs
#' in a \code{\link{phyloseq-class}} object.
#'
#' \code{UniFrac()} accesses the abundance
#' (\code{\link{otuTable-class}}) and a phylogenetic tree (\code{\link{phylo-class}})
#' data within an experiment-level (\code{\link{phyloseq-class}}) object.
#' If the tree and contingency table are separate objects, suggested solution
#' is to combine them into an experiment-level class
#' using the \code{\link{phyloseq}} function. For example, the following code
#'
#' \code{phyloseq(myOTUtable, myTree)}
#'
#' returns a \code{phyloseq}-class object that has been pruned and comprises
#' the minimum arguments necessary for \code{UniFrac()}. 
#'
#' Parallelization is possible for UniFrac calculated with the \code{\link{phyloseq-package}},
#' and is encouraged in the instances of large trees, many samples, or both.
#' Parallelization has been implemented via the \code{\link{foreach-package}}.
#' This means that parallel calls need to be preceded by 2 or more commands
#' that register the parallel ``backend''. This is acheived via your choice of
#' helper packages. One of the simplest seems to be the \emph{doParallel} package.
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
#' Furthermore, as of \code{R} version \code{2.14.0} and higher, a parallel package
#' is included as part of the core installation, \code{\link{parallel-package}}, 
#' and this can be used as the parallel backend with the \code{\link{foreach-package}}
#' using the adaptor package ``doParallel''. 
#' \url{http://cran.r-project.org/web/packages/doParallel/index.html}
#'
#' See the vignette for some simple examples for using doParallel.
#' \url{http://cran.r-project.org/web/packages/doParallel/vignettes/gettingstartedParallel.pdf}
#'
#' UniFrac-specific examples for doParallel are provided in the example
#' code below.  
#'
#' @usage UniFrac(physeq, weighted=FALSE, normalized=TRUE, parallel=FALSE, fast=TRUE)
#'
#' @param physeq (Required). \code{\link{phyloseq-class}}, containing at minimum
#'  a phylogenetic tree (\code{\link{phylo-class}}) and 
#'  contingency table (\code{\link{otuTable-class}}). See
#'  examples below for coercions that might be necessary.
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
#' @param fast (Optional). Logical. Do you want to use the ``Fast UniFrac''
#'  algorithm? Implemented natively in the \code{phyloseq-package}.
#'  This is the default and the recommended option. There should be no difference
#'  in the output between the two algorithms.
#'  Moreover, the original UniFrac algorithm
#'  only outperforms this implementation of fast-UniFrac if the datasets are so
#'  small 
#'  (approximated by the value of \code{nspecies(physeq) * nsamples(physeq)}) 
#'  that the difference in time is inconsequential (less than 1 second). In practice
#'  it does not appear that this parameter should ever be set to \code{FALSE}, but
#'  the option is nevertheless included in the package for comparisons and 
#'  instructional purposes.
#'
#' @return a sample-by-sample distance matrix, suitable for NMDS, etc.
#' 
#' @seealso \code{\link[vegan]{vegdist}}, \code{\link[picante]{unifrac}}
#'
#' @references \url{http://bmf.colorado.edu/unifrac/}
#'
#' The main implementation (Fast UniFrac) is adapted from the algorithm's
#' description in:
#' 
#' Hamady, Lozupone, and Knight,
#' ``Fast UniFrac: facilitating high-throughput phylogenetic analyses of 
#' microbial communities including analysis of pyrosequencing and PhyloChip data.''
#' The ISME Journal (2010) 4, 17--27.
#' 
#' \url{http://www.nature.com/ismej/journal/v4/n1/full/ismej200997a.html}
#'
#' See also additional descriptions of UniFrac in the following articles:
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
#' @examples
#' # ################################################################################
#' # # Perform UniFrac on esophagus data
#' # ################################################################################
#' # data("esophagus")
#' # (y <- UniFrac(esophagus, TRUE))
#' # UniFrac(esophagus, TRUE, FALSE)
#' # UniFrac(esophagus, FALSE)
#' # picante::unifrac(as(t(otuTable(esophagus)), "matrix"), tre(esophagus))
#' # ################################################################################
#' # # Try phylocom example data from picante package
#' # # It comes as a list, so you must construct the phyloseq object first.
#' # ################################################################################
#' # data("phylocom")
#' # (x1 <- phyloseq(otuTable(phylocom$sample, FALSE), phylocom$phylo))
#' # UniFrac(x1, TRUE)
#' # UniFrac(x1, TRUE, FALSE)
#' # UniFrac(x1, FALSE)
#' # picante::unifrac(phylocom$sample, phylocom$phylo)
#' # ################################################################################
#' # # Now try a parallel implementation using doParallel, which leverages the 
#' # # new 'parallel' core package in R 2.14.0+
#' # # Note that simply loading the 'doParallel' package is not enough, you must
#' # # call a function that registers the backend. In general, this is pretty easy
#' # # with the 'doParallel package' (or one of the alternative 'do*' packages)
#' # #
#' # # Also note that the esophagus example has only 3 samples, and a relatively small
#' # # tree. This is fast to calculate even sequentially and does not warrant
#' # # parallelized computation, but provides a good quick example for using UniFrac()
#' # # in a parallel fashion. The number of cores you should specify during the
#' # # backend registration, using registerDoParallel(), depends on your system and
#' # # needs. 3 is chosen here for convenience. If your system has only 2 cores, this
#' # # will probably fault or run slower than necessary.
#' # ################################################################################
#' # library(doParallel)
#' # data(esophagus)
#' # # For SNOW-like functionality (works on Windows):
#' # cl <- makeCluster(3)
#' # registerDoParallel(cl)
#' # UniFrac(esophagus, TRUE)
#' # # Force to sequential backed:
#' # registerDoSEQ()
#' # # For multicore-like functionality (will probably not work on windows),
#' # # register the backend like this:
#' # registerDoParallel(cores=3)
#' # UniFrac(esophagus, TRUE)
#' ################################################################################
setGeneric("UniFrac", function(physeq, weighted=FALSE, normalized=TRUE, parallel=FALSE, fast=TRUE){
	standardGeneric("UniFrac")
})
################################################################################
#' @aliases UniFrac,phyloseq-method
#' @rdname UniFrac-methods
setMethod("UniFrac", "phyloseq", function(physeq, weighted=FALSE, normalized=TRUE, parallel=FALSE, fast=TRUE){

	if( fast ){
		fastUniFrac(physeq, weighted, normalized, parallel)
	} else {
		originalUniFrac(physeq, weighted, normalized, parallel)
	}
	
})
################################################################################
################################################################################
# Fast UniFrac for R.
# Adapted from The ISME Journal (2010) 4, 17-27; doi:10.1038/ismej.2009.97;
# http://www.nature.com/ismej/journal/v4/n1/full/ismej200997a.html
################################################################################
#' @keywords internal
fastUniFrac <- function(physeq, weighted=FALSE, normalized=TRUE, parallel=FALSE){
	# # # require(picante); require(ape); require(foreach)

	# Access the needed components. Note, will error if missing in physeq.
	OTU  <- otuTable(physeq)
	tree <- tre(physeq)

	# Some important checks.
    if( is.null(tree$edge.length) ) {
        stop("Tree has no branch lengths, cannot compute UniFrac")
    }
    if( !is.rooted(tree) ) {
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

	# Make sure OTU is in species-are-rows orientation
	if( !speciesAreRows(OTU) ){OTU <- t(OTU)}    
   	# Enforce that tree and otuTable indices are the same order, 
   	# by re-ordering OTU, if needed
	if( !all(rownames(OTU) == species.names(tree)) ){
		OTU <- OTU[species.names(tree), ]
	}

	########################################
	# Build the requisite matrices as defined 
	# in the Fast UniFrac paper.
	########################################
	## This only needs to happen once in a call to UniFrac.
	## Notice that A and B do not appear in this section.
	
	# Begin by building the edge array (edge-by-sample matrix)
	edges <- 1:nrow(tree$edge)
	edge_array <- matrix(0, nrow=length(edges), ncol=nsamples(OTU), 
		dimnames=list(edge_index=edges, sample.names=sample.names(OTU)))

	# loop over each edge in the tree. Parallel version.
  	edge_array_list <- foreach( edge = edges, .packages="phyloseq") %dopar% {		
		# get the associated node in order to call internal2tips(). node number, n:
		n <- tree$edge[edge, 2]
		# get subset of tips (species) that descend from this branch
		edgeDescendants <- phyloseq:::internal2tips.self(tree, n, return.names = TRUE)
		# sum the descendants of this branch for every sample, assign to edge array (list)
		return( apply(OTU[edgeDescendants, ], 2, sum) )
	}
	# Quick loop to organize tree-edge results into a matrix, called edge_array
	for( edge in edges ){
		edge_array[edge, ] <- edge_array_list[[edge]]
	}
	
	# If unweighted-UniFrac, coerce to a presence-absence contingency, occ
	if(!weighted){
		# For unweighted UniFrac, convert the edge_array to an occurrence (presence/absence binary) array
		edge_occ <- (edge_array > 0) - 0
	}

	if( weighted & normalized ){
		# This is only relevant to weighted-UniFrac.
		# For denominator in the normalized distance, we need the age of each tip.
		# Get the tip ages from their associated edges (node.age gives the age of edges, ironically)
		### Note: this picante:node.age function is a computational bottleneck.
		### It could be vectorized/parallelized 
		tipAges <- picante::node.age(tree)$ages[ which(tree$edge[, 2] %in% 1:length(tree$tip.label)) ]
		names(tipAges) <- tree$tip.label		
	}

	########################################	
   	# optionally-parallel implementation with foreach
	########################################
  	distlist <- foreach( i = spn, .packages="phyloseq") %dopar% {
		A  <- i[1]
		B  <- i[2]
		AT <- sampleSums(OTU)[A]
		BT <- sampleSums(OTU)[B]		
		if( weighted ){ # weighted UniFrac
			# subset matrix to just columns A and B
			edge_array_AB <- edge_array[, c(A, B)]
			# Perform UFwi equivalent, "inline" with apply on edge_array_AB
			wUF_branchweight <- apply(edge_array_AB, 1, function(br, A, B, ATBT){
				abs((br[A]/ATBT[A]) - (br[B]/ATBT[B]))
			}, A, B, c(AT, BT))
			# calculate the w-UF numerator
			numerator <- sum(tree$edge.length * wUF_branchweight)
			# if not-normalized weighted UniFrac, just return "numerator";
			# the u-value in the w-UniFrac description
			if(!normalized){
				return(numerator)
			} else {
				# denominator (assumes tree-indices and otuTable indices are same order)
				denominator <- sum( tipAges * (OTU[, A]/AT + OTU[, B]/BT) )
				# return the normalized weighted UniFrac values
				return(numerator / denominator)
			}
		} else { # unweighted UniFrac
			# subset matrix to just columns A and B
			edge_occ_AB <- edge_occ[, c(A, B)]
			# keep only the unique branches, sum the lengths
			edge_uni_AB_sum <- sum( (tree$edge.length * edge_occ_AB)[apply(edge_occ_AB, 1, sum) < 2, ] )
			# Normalize this sum to the total branches among these two samples, A and B
			uwUFpairdist <- edge_uni_AB_sum / sum( tree$edge.length[apply(edge_occ_AB, 1, sum) > 0] )
			return(uwUFpairdist)
		}
	}
	
    # This is in serial, but it is quick.
    distlist2distmat <- function(i, spn, DL){
    	UniFracMat[ spn[[i]][2], spn[[i]][1] ] <<- DL[[i]]
    }
	junk <- sapply(1:length(spn), distlist2distmat, spn, distlist)
    
    return(as.dist(UniFracMat))
	
}
################################################################################