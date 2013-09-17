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
# @aliases vegdist,otu_table-method
# @rdname vegdist-methods
setMethod("vegdist", "otu_table", function(x, method = "bray", binary = FALSE,
	diag = FALSE, upper = FALSE, na.rm = FALSE, ...){
	# Make sure in sample-by-species orientation
	if( taxa_are_rows(x) ){x <- t(x)}
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
	# Simply access the otu_table
	x <- otu_table(x)
	vegdist(x, method, binary, diag, upper, na.rm, ...)	
})
################################################################################
# @importFrom vegan estimateR
# @importFrom vegan diversity
################################################################################
#' Summarize alpha diversity
#'
#' Performs a number of standard alpha diversity estimates, 
#' and returns the results as a \code{data.frame}.
#' Strictly speaking, this function is not only estimating richness,
#' despite its name. 
#' It can operate on the cumulative population of all
#' samples in the dataset, or by repeating the richness estimates for each
#' sample individually.
#' NOTE: You must use untrimmed datasets
#' for meaningful results, as these estimates (and even the ``observed'' richness)
#' are highly dependent on the number of singletons. You can always trim the data
#' later on if needed, just not before using this function.
#' 
#' @param physeq (Required). \code{\link{phyloseq-class}}, or alternatively, 
#'  an \code{\link{otu_table-class}}. The data about which you want to estimate
#'  the richness.
#'
#' @param split (Optional). Logical. Should a separate set of richness estimates
#'  be performed for each sample? Or alternatively, pool all samples and 
#'  estimate richness of the entire set.
#'  
#' @param measures (Optional). Default is \code{NULL}, meaning that
#'  all available alpha-diversity measures will be included.
#'  Alternatively, you can specify one or more measures
#'  as a character vector of measure names.
#'  Values must be among those supported:
#'  \code{c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")}.
#'
#' @return A \code{data.frame} of the richness estimates, and their standard error.
#' 
#' @seealso 
#'  Check out the custom plotting function, \code{\link{plot_richness}},
#'  for easily showing the results of different estimates, 
#'  with method-specific error-bars.
#'  Also check out the internal functions borrowed from the \code{vegan} package:
#'  
#'  \code{\link[vegan]{estimateR}}
#'  
#'  \code{\link[vegan]{diversity}}
#'  
#'  \code{\link[vegan]{fisherfit}}
#'
#' @import vegan
#' @export
#' @examples 
#' ## There are many more interesting examples at the phyloseq online tutorials.
#' ## http://joey711.github.com/phyloseq/plot_richness-examples
#'  data("soilrep")
#'  # Default is all available measures
#'  estimate_richness(soilrep)
#'  # Specify just one:
#'  estimate_richness(soilrep, measures="Observed")
#'  # Specify a few:
#'  estimate_richness(soilrep, measures=c("Observed", "InvSimpson", "Shannon", "Chao1"))
estimate_richness <- function(physeq, split=TRUE, measures=NULL){

	if( !any(otu_table(physeq)==1) ){
	  # Check for singletons, and then warning if they are missing.
	  # These metrics only really meaningful if singletons are included.
		warning(
			"The data you have provided does not have\n",
			"any singletons. This is highly suspicious. Results of richness\n",
			"estimates (for example) are probably unreliable, or wrong, if you have already\n",
			"trimmed low-abundance taxa from the data.\n",
			"\n",
			"We recommended that you find the un-trimmed data and retry."
		)
	}
	
	# If we are not splitting sample-wise, sum the species. Else, enforce orientation.
	if( !split ){
		OTU <- taxa_sums(physeq)		
	} else if( split ){
		OTU <- as(otu_table(physeq), "matrix")
		if( taxa_are_rows(physeq) ){ OTU <- t(OTU) }
	}
	
	# Define renaming vector:
	renamevec = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")
	names(renamevec) <- c("S.obs", "S.chao1", "S.ACE", "shannon", "simpson", "invsimpson", "fisher")
	# If measures was not explicitly provided (is NULL), set to all supported methods
	if( is.null(measures) ){
	  measures = as.character(renamevec)
	}
  # Rename measures if they are in the old-style
  if( any(measures %in% names(renamevec)) ){
    measures[measures %in% names(renamevec)] <- renamevec[names(renamevec) %in% measures]
  }
  
  # Stop with error if no measures are supported
  if( !any(measures %in% renamevec) ){
    stop("None of the `measures` you provided are supported. Try default `NULL` instead.")
  }
  
  # Initialize to NULL
  outlist = vector("list")
	# Some standard diversity indices
  estimRmeas = c("Chao1", "Observed", "ACE")
	if( any(estimRmeas %in% measures) ){ 
    outlist <- c(outlist, list(t(data.frame(estimateR(OTU)))))
	}
	if( "Shannon" %in% measures ){
    outlist <- c(outlist, list(shannon = diversity(OTU, index="shannon")))
	}
	if( "Simpson" %in% measures ){
	  outlist <- c(outlist, list(simpson = diversity(OTU, index="simpson")))
	}
	if( "InvSimpson" %in% measures ){
	  outlist <- c(outlist, list(invsimpson = diversity(OTU, index="invsimpson")))
	}
	if( "Fisher" %in% measures ){
    fisher = tryCatch(fisher.alpha(OTU, se=TRUE), 
      warning=function(w){
        warning("phyloseq::estimate_richness: Warning in fisher.alpha(). See `?fisher.fit` or ?`fisher.alpha`. Treat fisher results with caution")
        suppressWarnings(fisher.alpha(OTU, se=TRUE)[, c("alpha", "se")])
      }
    )
	  colnames(fisher) <- c("Fisher", "se.fisher")
	  outlist <- c(outlist, list(fisher))
	}
  out = do.call("cbind", outlist)
  # Rename columns per renamevec
  namechange = intersect(colnames(out), names(renamevec))
  colnames(out)[colnames(out) %in% namechange] <- renamevec[namechange]
  # Final prune to just those columns related to "measures". Use grep.
  colkeep = sapply(paste0("(se\\.){0,}", measures), grep, colnames(out), ignore.case=TRUE)
  out = out[, sort(unique(unlist(colkeep))), drop=FALSE]
  # Make sure that you return a data.frame for reliable performance.
  out <- as.data.frame(out)
	return(out)
}
################################################################################
