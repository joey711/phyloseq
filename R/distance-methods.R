################################################################################
#' General distance / dissimilarity index calculator
#'
#' Takes a \code{\link{phyloseq-class}} object and method option, and returns
#'  a \code{\link{dist}}ance object suitable for certain 
#'  ordination methods and other distance-based analyses. 
#'  There are currently 44 explicitly supported method options, as well as
#'  user-provided arbitrary methods via an interface to
#'  \code{\link{designdist}}. For the complete list of currently
#'  supported options/arguments to the \code{method} parameter, 
#'  type \code{distance("list")} at the command-line.
#'  Only
#'  sample-wise distances are currently supported (the \code{type} argument),
#'  but eventually species-wise (OTU-wise)
#'  distances will be supported as well. 
#'
#'  Depending on the \code{method}
#'  argument, \code{distance()} wraps one of 
#'  \code{\link{UniFrac}},
#'  \code{\link{DPCoA}},
#'  \code{\link{JSD}},
#'  \code{\link[vegan]{vegdist}},
#'  \code{\link[vegan]{betadiver}},
#'  \code{\link[vegan]{designdist}}, or
#'  \code{\link{dist}}.
#' 
#' @param physeq (Required).  A \code{\link{phyloseq-class}} or
#'  an \code{\link{otu_table-class}} object. The latter is only appropriate
#'  for methods that do not require any additional data (one-table). 
#'  For example, the ``unifrac'' option (\code{\link{UniFrac}}) requires
#'  \code{\link{phyloseq-class}} that contains both an \code{otu_table}
#'  and a phylogenetic tree (\code{phylo}).
#'
#' @param method (Optional). A character string. Default is \code{"unifrac"}.
#'  Provide one of the 44 currently supported options. 
#'  To see a list of supported options, enter the following into the command line:
#' 
#'  \code{distance("list")}
#'  
#'  For further details and additional arguments,
#'  see the documentation for the supprting functions, linked below
#'  under ``See Also''. 
#'  
#'  In particular, there are three methods included
#'  by the \code{\link{phyloseq-package}}, and accessed by the following
#'  \code{method} options:
#' 
#'  \code{"unifrac"}, for UniFrac based distances, \code{\link{UniFrac}};
#'
#'  \code{"dpcoa"}, sample-wise distance from Double Principle 
#'   Coordinate Analysis, \code{\link{DPCoA}};
#' 
#'  \code{"jsd"}, for Jensen-Shannon Divergence, \code{\link{JSD}}; 
#'  
#'  and it is recommended that you see their documentation
#'  for details, references, background and examples for use. 
#'
#'  Alternatively, you can provide
#'  a character string that defines a custom distance method, if it has the form
#'  described in \code{\link{designdist}}. 
#'
#' @param type (Optional). A character string. The type of pairwise comparisons
#'  being calculated: sample-wise or taxa-wise. The default is 
#'  \code{c("samples")}.
#'
#' @param ... Additional arguments passed on to the appropriate distance 
#'  function, determined by the \code{method} argument.
#'
#' @return An object of class ``\code{\link{dist}}'' suitable for certain 
#'  ordination methods and other distance-based analyses.
#' 
#' @seealso 
#'  \code{\link{plot_ordination}}, 
#'  \code{\link{UniFrac}},
#'  \code{\link{DPCoA}},
#'  \code{\link{JSD}},
#'  \code{\link[vegan]{vegdist}},
#'  \code{\link[vegan]{betadiver}},
#'  \code{\link[vegan]{designdist}},
#'  \code{\link{dist}}.
#'
#' @importFrom vegan betadiver
#' @importFrom vegan designdist
#' @importFrom vegan vegdist
#' @export
#' @examples 
#' data(esophagus)
#' distance(esophagus) # Unweighted UniFrac
#' distance(esophagus, weighted=TRUE) # weighted UniFrac
#' distance(esophagus, "jaccard") # vegdist jaccard
#' distance(esophagus, "gower") # vegdist option "gower"
#' distance(esophagus, "g") # designdist method option "g"
#' distance(esophagus, "minkowski") # invokes a method from the base dist() function.
#' distance(esophagus, "(A+B-2*J)/(A+B)") # designdist custom distance
#' distance("help")
#' distance("list")
#' help("distance")
distance <- function(physeq, method="unifrac", type="samples", ...){
	# # Can't do partial matching directly, because too many similar options.
	# # Determine if method argument matches any options exactly.
	# # If not, call designdist
	vegdist_methods <- c("manhattan", "euclidean", "canberra", "bray", 
		"kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", 
		"mountford", "raup" , "binomial", "chao", "cao")

	# Special methods (re)defined in phyloseq
	phyloseq_methods <- c("unifrac", "dpcoa", "jsd")
	
	# The standard distance methods
	dist_methods <- c("maximum", "binary", "minkowski")
	# Only keep the ones that are NOT already in vegdist_methods
	dist_methods <- dist_methods[!dist_methods %in% intersect(vegdist_methods, dist_methods)]

	# The methods supported by vegan::betadiver function.
	betadiver_methods <- c("w", "-1", "c", "wb", "r", "I", "e", "t", "me", "j",
		"sor", "m", "-2", "co", "cc", "g", "-3", "l", "19", "hk", "rlb",
		"sim", "gl", "z")
	
	method.list <- list(
		UniFrac    = "unifrac",
		DPCoA      = "dpcoa",
		JSD        = "jsd",
		vegdist    = vegdist_methods,
		betadiver  = betadiver_methods,
		dist       = dist_methods,
		designdist = "ANY"
	)
	
	if(class(physeq) == "character"){
		if( physeq=="help" ){
			cat("Available arguments to methods:\n")
			print(method.list)
			cat("Please be exact, partial-matching not supported.\n")
			cat("Can alternatively provide a custom distance.\n")
			cat("See:\n help(\"distance\") \n")
			return()
		}
		if( physeq=="list" ){
			return(c(method.list))
		}		
	}
  
	# Define the function call to build
	if( method %in% phyloseq_methods ){
		if( method == "unifrac"){ return(UniFrac(physeq, ...)) }
		if( method == "jsd"    ){ return(    JSD(physeq, ...)) }
		if( method == "dpcoa"  ){ return(dist(DPCoA(physeq, ...)$RaoDis)) }
	} else if( method %in% vegdist_methods ){
		dfun <- "vegdist"
	} else if( method %in% betadiver_methods ){
		dfun <- "betadiver"
	} else if( method %in% dist_methods ){
		dfun <- "dist"		
	} else {
		dfun <- "designdist"
	}
	
	# get the extra arguments to pass to functions (this can be empty)
	extrargs <- list(...)	

	# # non-phyloseq methods are assumed to be based on otu_table-only (for now)
	# # If necessary (non phyloseq funs), enforce orientation, build function.
	OTU <- otu_table(physeq)

	# disambiguate type argument... Must be "species" for vegan integration...
	if( type %in% c("taxa", "species", "OTUs", "otus", "otu") ){
		type <- "species"
	}

	# Test type, and enforce orientation accordingly	
	if( type == "species"){
		# For species-distance, species need to be rows (vegan-style)
		if( !taxa_are_rows(OTU) ){OTU <- t(OTU)}	
	} else if( type == "samples" ){
		# For sample-distance, samples need to be rows (vegan-style)
		if( taxa_are_rows(OTU) ){OTU <- t(OTU)}
	} else {
		stop("type argument must be one of \n (1) samples \n or \n (2) species")
	}	

	OTU <- as(OTU, "matrix")
	fun.args <- c(list(OTU, method=method), extrargs)	
	return( do.call(dfun, fun.args) )
} 
################################################################################
################################################################################
# Shannon-Jensen Divergence, in R.
################################################################################
#' @keywords internal
JSD.pair <- function(x, y){
	###Function to compute Shannon-Jensen Divergence
	###x and y are the frequencies for the same p categories
	u <- x/sum(x)
	v <- y/sum(y)
	m <- (u+v)/2
	if (all(u*v>0)){
		d <- (u*log(u/m)+v*log(v/m))/2
	} else {
		P1 <- u*log(u/m)
		P2 <- v*log(v/m)
		P1[is.nan(P1)] <- 0
		P2[is.nan(P2)] <- 0
		d <- (P1+P2)/2
	}
	return(sum(d))
}
################################################################################
#' Calculate the Jensen-Shannon Divergence (distance)
#'
#' This is a phyloseq-specific implementation of the Jensen-Shannon Divergence
#' for comparing pairs of microbial communities (samples) in an experiment. 
#' The expectation is that you have many samples (say. more than two) and you
#' want a distance matrix on which will perform further analysis. \code{JSD} is
#' intended to be ``wrapped'' by the more general \code{\link{distance}} 
#' function in phyloseq, and it can be invoked using \code{"jsd"} as the 
#' argument to the \code{method} parameter of \code{\link{distance}}.
#'
#' One of the motivations for providing JSD in phyloseq was its recent use in 
#' the analysis of the \code{\link{enterotype}} dataset. 
#'
#' @usage JSD(physeq, parallel=FALSE)
#' 
#' @param physeq (Required). \code{\link{phyloseq-class}}. 
#'  The phyloseq data on which to compute the 
#'  pairwise sample distance matrix.
#'
#' @param parallel (Optional). Logical. Default \code{FALSE}. Should the 
#'  calculation be run using parallel processors? Note that this only actually
#'  runs in parallel if you have registered a parallel ``backend'' prior to 
#'  calling this function. Further details are provided in the documentation
#'  for \code{\link{UniFrac}}.
#'
#' @return An object of class ``\code{\link{dist}}'' suitable for certain 
#'  ordination methods and other distance-based analyses.
#'  See \code{\link{distance}}.
#'  
#' @seealso
#'  \code{\link{distance}} 
#' 
#'  \code{\link{enterotype}}
#' 
#'  \url{http://en.wikipedia.org/wiki/Jensen-Shannon_divergence}
#'
#' @references
#' Jensen-Shannon Divergence and Hilbert space embedding.  
#' Bent Fuglede and Flemming Topsoe University of Copenhagen,
#' Department of Mathematics
#' \url{http://www.math.ku.dk/~topsoe/ISIT2004JSD.pdf}
#'
#' @author
#'  Susan Holmes \email{susan@@stat.stanford.edu}.
#'  Adapted for phyloseq by Paul J. McMurdie.
#'
#' @keywords internal
#' @examples
#' # library(doParallel)  # Do this and next line only if you have multi-cores
#' # registerDoParallel(cores=6)
#' # data(enterotype)
#' # # ent.jsd <- JSD(enterotype, TRUE) # internal only
#' # ent.jsd <- distance(enterotype, "jsd", parallel=TRUE)
#' # ent.PCoA <- ordinate(enterotype, "PCoA", ent.jsd) # Perform principle coordinate analysis
#' # p <- plot_ordination(enterotype, ent.PCoA, color="Enterotype", shape="SeqTech") 
#' # (p <- p + geom_point(size=5, alpha=0.5))
JSD <- function(physeq, parallel=FALSE){
	OTU <- otu_table(physeq)
	### Some parallel-foreach housekeeping.    
	# If user specifies not-parallel run (the default), register the sequential "back-end"
	if( !parallel ){ registerDoSEQ() }
	# create N x 2 matrix of all pairwise combinations of samples.
	spn <- combn(sample_names(physeq), 2, simplify=FALSE)
	# initialize DistMat with NAs
	DistMat <- matrix(NA, nsamples(physeq), nsamples(physeq))
	# define the rows/cols of DistMat with the sample names (rownames)    
	rownames(DistMat) <- sample_names(physeq)
	colnames(DistMat) <- sample_names(physeq)
	## Format coercion
	# Coerce to the vegan orientation, with species as columns
	if(taxa_are_rows(physeq)){ OTU <- t(OTU) }
	# Coerce OTU to matrix for calculations.
	OTU <- as(OTU, "matrix")
	# optionally-parallel implementation with foreach
	distlist <- foreach( i = spn, .packages="phyloseq") %dopar% {
	  A <- i[1]
	  B <- i[2]
	  return( phyloseq:::JSD.pair(OTU[A, ], OTU[B, ]) )
	}
	# return(distlist)
	# This is in serial, but it is quick.
	distlist2distmat <- function(i, spn, DL){
	  DistMat[ spn[[i]][2], spn[[i]][1] ] <<- DL[[i]]
	}
	junk <- sapply(1:length(spn), distlist2distmat, spn, distlist)	
	return(as.dist(DistMat))
}
##############################################################################
#' Calculate weighted or unweighted (Fast) UniFrac distance for all sample pairs.
#'
#' This function calculates the (Fast) UniFrac distance for all sample-pairs
#' in a \code{\link{phyloseq-class}} object.
#'
#' \code{UniFrac()} accesses the abundance
#' (\code{\link{otu_table-class}}) and a phylogenetic tree (\code{\link{phylo-class}})
#' data within an experiment-level (\code{\link{phyloseq-class}}) object.
#' If the tree and contingency table are separate objects, suggested solution
#' is to combine them into an experiment-level class
#' using the \code{\link{phyloseq}} function. For example, the following code
#'
#' \code{phyloseq(myotu_table, myTree)}
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
#'  contingency table (\code{\link{otu_table-class}}). See
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
#' @param fast (Optional). Logical. DEPRECATED.
#'  Do you want to use the ``Fast UniFrac''
#'  algorithm? Implemented natively in the \code{phyloseq-package}.
#'  \code{TRUE} is now the only supported option.
#'  There should be no difference in the output between the two algorithms.
#'  Moreover, the original UniFrac algorithm
#'  only outperforms this implementation of fast-UniFrac if the datasets are so
#'  small 
#'  (approximated by the value of \code{ntaxa(physeq) * nsamples(physeq)}) 
#'  that the difference in time is inconsequential (less than 1 second).
#'  In practice it does not appear that this parameter should 
#'  have ever been set to \code{FALSE}, and therefore 
#'  the original UniFrac implementation perhaps never should have been supported here.
#'  For legacy code support the option is now deprecated here
#'  (the implementation was an internal function, anyway)
#'  and the \code{fast} option will remain for one release cycle before
#'  being removed completely.
#'
#' @return a sample-by-sample distance matrix, suitable for NMDS, etc.
#' 
#' @seealso
#' 
#' \code{\link{distance}}
#' 
#' \code{unifrac} in the picante package.
#'
#' @references
#' 
#' \url{http://bmf.colorado.edu/unifrac/}
#'
#' The main implementation (Fast UniFrac) is adapted from the algorithm's
#' description in:
#' 
#' Hamady, Lozupone, and Knight,
#' ``\href{http://www.nature.com/ismej/journal/v4/n1/full/ismej200997a.html}{Fast UniFrac:}
#' facilitating high-throughput phylogenetic analyses of 
#' microbial communities including analysis of pyrosequencing and PhyloChip data.''
#' The ISME Journal (2010) 4, 17--27.
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
#' ################################################################################
#' # Perform UniFrac on esophagus data
#' ################################################################################
#' data("esophagus")
#' (y <- UniFrac(esophagus, TRUE))
#' UniFrac(esophagus, TRUE, FALSE)
#' UniFrac(esophagus, FALSE)
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
#' @importFrom ape is.rooted
#' @importFrom ape root
setMethod("UniFrac", "phyloseq", function(physeq, weighted=FALSE, normalized=TRUE, parallel=FALSE, fast=TRUE){
  if(is.null(phy_tree(physeq)$edge.length)){
    stop("Tree has no branch lengths. See tree$edge.length. Cannot compute UniFrac without branch lengths")
  }  
  # Check if tree is rooted, set random root with warning if it is not.
  if( !is.rooted(phy_tree(physeq)) ){
    #phy_tree(physeq) <- unroot(phy_tree(physeq))
    randoroot = sample(taxa_names(physeq), 1)
    warning("Randomly assigning root as -- ", randoroot, " -- in the phylogenetic tree in the data you provided.")
    phy_tree(physeq) <- root(phy=phy_tree(physeq), outgroup=randoroot, resolve.root=TRUE, interactive=FALSE)
    if( !is.rooted(phy_tree(physeq)) ){
      stop("Problem automatically rooting tree. Make sure your tree is rooted before attempting UniFrac calculation. See ?ape::root")
    }
  } 
	if( fast ){
		fastUniFrac(physeq, weighted, normalized, parallel)
	} else {
    warning("Option `fast=FALSE` is deprecated. Only 'fast' UniFrac is supported in phyloseq.")
	  fastUniFrac(physeq, weighted, normalized, parallel)
	}
})
################################################################################
# Fast UniFrac for R.
# Adapted from The ISME Journal (2010) 4, 17-27; doi:10.1038/ismej.2009.97;
# http://www.nature.com/ismej/journal/v4/n1/full/ismej200997a.html
################################################################################
#' @importFrom ape prop.part
#' @importFrom ape reorder.phylo
#' @keywords internal
#' @import foreach
fastUniFrac <- function(physeq, weighted=FALSE, normalized=TRUE, parallel=FALSE){
	# Access the needed components. Note, will error if missing in physeq.
	OTU  <- otu_table(physeq)
	tree <- phy_tree(physeq)
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
	spn <- combn(sample_names(physeq), 2, simplify=FALSE)
	# Make sure OTU is in species-are-rows orientation
	if( !taxa_are_rows(physeq) ){OTU <- t(OTU)}
  # Convert to standard matrix
	OTU <- as(OTU, "matrix")  
	# Enforce that tree and otu_table indices are the same order, 
	# by re-ordering OTU, if needed
	if( !all(rownames(OTU) == taxa_names(tree)) ){
	  OTU <- OTU[taxa_names(tree), ]
	}
	########################################
	# Build the requisite matrices as defined 
	# in the Fast UniFrac article.
	########################################
	## This only needs to happen once in a call to UniFrac.
	## Notice that A and B do not appear in this section.
	# Begin by building the edge descendants matrix (edge-by-sample)
  # `edge_array`
  #
	# Create a list of descendants, starting from the first internal node (root)
	descList <- prop.part(tree, check.labels = FALSE)
	# Add the terminal edge descendants (tips). By definition, can only have one descendant
	descList <- c(as.list(1:length(tree$tip.label)), descList)
	# Convert `descList` to `edge_array` that matches the order of things in `tree$edge`
	edge_array <- matrix(0, nrow=nrow(tree$edge), ncol=nsamples(physeq), 
	                      dimnames=list(NULL, sample_names=sample_names(physeq)))
	for(i in 1:nrow(tree$edge)){
    # For each entry in the tree$edge table, sum the descendants for each sample
    # `tree$edge[i, 2]` is the node ID.
	  edge_array[i, ] <- colSums(OTU[descList[[tree$edge[i, 2]]], , drop=FALSE], na.rm = TRUE)
	}
  # Remove unneeded variables. `descList` in particular could be large-ish.
  rm(descList)
	# If unweighted-UniFrac, coerce to a presence-absence contingency, occ
	if(!weighted){
		# For unweighted UniFrac, convert the edge_array to an occurrence (presence/absence binary) array
		edge_occ <- (edge_array > 0) - 0
	}
	if( weighted & normalized ){
		# This is only relevant to weighted-UniFrac.
		# For denominator in the normalized distance, we need the age of each tip.
	  # 'z' is the tree in postorder order used in calls to .C
	  # Descending order of left-hand side of edge (the ancestor to the node)
	  z = reorder.phylo(tree, order="postorder")
	  # Call phyloseq-internal function that in-turn calls ape's internal
	  # horizontal position function, in C, using the re-ordered phylo object, `z`
	  tipAges = ape_node_depth_edge_length(Ntip = length(tree$tip.label),
	                                       Nnode = tree$Nnode, 
	                                       edge = z$edge, 
	                                       Nedge = nrow(tree$edge)[1],
	                                       edge.length = z$edge.length)
	  # Keep only the tips, and add the tip labels in case `z` order differs from `tree`
	  tipAges <- tipAges[1:length(tree$tip.label)]
	  names(tipAges) <- z$tip.label	
    # Explicitly re-order tipAges to match OTU
	  tipAges <- tipAges[rownames(OTU)]
	}
	########################################	
  # optionally-parallel implementation with foreach
	########################################
	samplesums = sample_sums(physeq)
	distlist <- foreach( i = spn, .packages="phyloseq") %dopar% {
	  A  <- i[1]
	  B  <- i[2]
	  AT <- samplesums[A]
	  BT <- samplesums[B]
	  if( weighted ){
      # weighted UniFrac
	    wUF_branchweight <- abs(edge_array[, A]/AT - edge_array[, B]/BT)
	    # calculate the w-UF numerator
	    numerator <- sum({tree$edge.length * wUF_branchweight}, na.rm = TRUE)
	    # if not-normalized weighted UniFrac, just return "numerator";
	    # the u-value in the w-UniFrac description
	    if(!normalized){
	      return(numerator)
	    } else {
	      # denominator (assumes tree-indices and otu_table indices are same order)
	      denominator <- sum({tipAges * (OTU[, A]/AT + OTU[, B]/BT)}, na.rm = TRUE)
	      # return the normalized weighted UniFrac values
	      return(numerator / denominator)
	    }
	  } else {
      # Unweighted UniFrac
	    # Subset matrix to just columns A and B
	    edge_occ_AB <- edge_occ[, c(A, B)]
      # Keep only the unique branches. Sum the lengths
      edge_uni_AB_sum <- sum((tree$edge.length * edge_occ_AB)[rowSums(edge_occ_AB, na.rm=TRUE) < 2, ], na.rm=TRUE)
	    # Normalize this sum to the total branches among these two samples, A and B
	    uwUFpairdist <- edge_uni_AB_sum / sum(tree$edge.length[rowSums(edge_occ_AB, na.rm=TRUE) > 0])
	    return(uwUFpairdist)
	  }
	}
	# Initialize UniFracMat with NAs
	UniFracMat <- matrix(NA_real_, nsamples(physeq), nsamples(physeq))
	rownames(UniFracMat) <- colnames(UniFracMat) <- sample_names(physeq)
  # Matrix-assign lower-triangle of UniFracMat. Then coerce to dist and return.
  	matIndices <- do.call(rbind, spn)[, 2:1]
  	# Take care of edge case where there are two samples -> 1 pair of indices -> rbind doesn't return a matrix
  	if(!is.matrix(matIndices)) matIndices <- matrix(matIndices, ncol=2)
	UniFracMat[matIndices] <- unlist(distlist)
	return(as.dist(UniFracMat))	
}
################################################################################
