################################################################################
# Define S3 methods for scores (originally defined by vegan-package) 
# to work for other ordination results
# vegan:::scores.default
################################################################################
# pcoa-class, from pcoa{ape}
#' @import ape
#' @keywords internal
scores.pcoa <- function(x, choices=NULL, display="sites", ...){
	if(is.null(choices)){
		choices <- colnames(x$vectors)
	}
	if( !display %in% c("sites", "samples") ){
		# MDS/PCoA only provides coordinates of the elements in the
		# distance matrix, usually sites/samples, so species (etc.)
		# not an option...
		return( NULL )
	} else {
		return( x$vectors[, choices] )		
	}
}
# dpcoa-class, from ade4::dpcoa or phyloseq::DPCoA
# @import ade4
#' @keywords internal
scores.dpcoa <- function(x, choices=NULL, display="sites", ...){
	ifelse(display=="species", coords <- x$l1, coords <- x$l2)
	if( is.null(choices) ){
		choices <- colnames(coords)
	}
	return( coords[, choices] )
}
################################################################################
################################################################################
# Extend vegdist for phyloseq classes
# @importFrom vegan vegdist
################################################################################
# \code{\link[vegan]{vegdist}} wrapper for phyloseq classes
#
# Trivially-extended S4 method from the \code{\link[vegan]{vegdist}} function,
# such that S4 classes from the \code{\link{phyloseq-package}} are properly
# handled / accessed. All parameters passed on to \code{\link[vegan]{vegdist}}
# verbatim.
#
# @seealso \code{\link[vegan]{vegdist}} 
# @import vegan
# @rdname vegdist-methods
# @docType methods
# @aliases vegdist
#
# @examples
# # data(esophagus)
# # vegdist(esophagus, "jaccard")
#' @keywords internal
setGeneric("vegdist")
################################################################################
# @aliases vegdist,otuTable-method
# @rdname vegdist-methods
setMethod("vegdist", "otuTable", function(x, method = "bray", binary = FALSE,
	diag = FALSE, upper = FALSE, na.rm = FALSE, ...){
	# Make sure in sample-by-species orientation
	if( speciesAreRows(x) ){x <- t(x)}
	# Convert to simple matrix
	x <- as(x, "matrix")
	# pass to standard method (compiled C)
	vegdist(x, method, binary, diag, upper, na.rm, ...)	
})
################################################################################
# @aliases vegdist,phyloseq-method
# @rdname vegdist-methods
setMethod("vegdist", "phyloseq", function(x, method = "bray", binary = FALSE,
	diag = FALSE, upper = FALSE, na.rm = FALSE, ...){
	# Simply access the otuTable
	x <- otuTable(x)
	vegdist(x, method, binary, diag, upper, na.rm, ...)	
})
################################################################################
# @importFrom vegan estimateR
# @importFrom vegan diversity
################################################################################
#' Summarize richness estimates
#'
#' Performs a number of standard richness estimates, and returns the results
#' as a \code{data.frame}. Can operate on the cumulative population of all
#' samples in the dataset, or by repeating the richness estimates for each
#' sample individually.
#' NOTE: You must use untrimmed datasets
#' for meaningful results, as these estimates (and even the ``observed'' richness)
#' are highly dependent on the number of singletons. You can always trim the data
#' later on if needed, just not before using this function.
#'
#' @usage estimate_richness(physeq, split=TRUE)
#' 
#' @param physeq (Required). \code{\link{phyloseq-class}}, or alternatively, 
#'  an \code{\link{otuTable-class}}. The data about which you want to estimate
#'  the richness.
#'
#' @param split (Optional). Logical. Should a separate set of richness estimates
#'  be performed for each sample? Or alternatively, pool all samples and 
#'  estimate richness of the entire set.
#'
#' @return A \code{data.frame} of the richness estimates, and their standard error.
#' 
#' @seealso 
#'  Check out the custom plotting function, \code{\link{plot_richness_estimates}},
#'  for easily showing the results of different estimates, with method-specific
#'  error-bars. Also check out the internal functions borrowed from the \code{vegan}
#'  package:
#'  \code{\link[vegan]{estimateR}},
#'  \code{\link[vegan]{diversity}}
#'
#' @import vegan
#' @export
#' @examples 
#'  data(GlobalPatterns)
#'  ( S.GP <- estimate_richness(GlobalPatterns) )
#'  # # Make the plots
#'  # plot_richness_estimates(GlobalPatterns, "SampleType")
#'  # plot_richness_estimates(GlobalPatterns, "SampleType", "SampleType")
#'  # For more plotting examples, see plot_richness_estimates()
#' 
estimate_richness <- function(physeq, split=TRUE){
	# Check for singletons, and then warning if they are missing.
	# These metrics only really meaningful if singletons are included.
	if( !any(otuTable(physeq)==1) ){
		warning("The experiment object you have provided does not have\n",
		"any singletons. This is highly suspicious. Results of richness\n",
		"estimates are probably unreliable, or wrong, if you have already\n",
		"trimmed low-abundance taxa from the data.\n",
		"\n",
		"It is recommended that you find the un-trimmed data and retry.",
		)
	}
	
	# If we are not splitting sample-wise, sum the species. Else, enforce orientation.
	if( !split ){
		OTU <- speciesSums(physeq)		
	} else if( split ){
		OTU <- as(otuTable(physeq), "matrix")
		if( speciesAreRows(physeq) ){ OTU <- t(OTU) }
	}
	
	# Some standard richness parameters
	richness <- round(estimateR(OTU))
	shannon	 <- round(diversity( OTU ), 2)
	simpson  <- round(diversity( OTU, index="simpson"), 2)
	# # fisher   <- round(fisher.alpha( OTU))
	
	return( t(rbind(richness, shannon, simpson)) )
}
################################################################################