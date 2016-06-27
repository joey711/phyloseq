################################################################################
#' Perform an ordination on phyloseq data
#'
#' This function wraps several commonly-used ordination methods. The type of 
#' ordination depends upon the argument to \code{method}. Try 
#' \code{ordinate("help")} or \code{ordinate("list")} for the currently
#' supported method options.
#'
#' @param physeq (Required). Phylogenetic sequencing data
#'  (\code{\link{phyloseq-class}}). The data on which you want to perform
#'  the ordination. In general, these methods will be based in some fashion on
#'  the abundance table ultimately stored as a contingency matrix 
#'  (\code{\link{otu_table-class}}). If you're able to import data into 
#'  \code{\link{phyloseq-class}} format, than you don't need to worry, as an
#'  \code{otu_table} is a required component of this class. In addition, some 
#'  ordination methods require additional data, like a constraining variable
#'  or phylogenetic tree. If that is the case, the relevant data should be
#'  included in \code{physeq} prior to running. Integrating the data in this way
#'  also results in these different data components being checked for validity
#'  and completeness by the method.
#'
#' @param method (Optional). A character string. Default is \code{"DCA"}.
#'
#'  Currently supported method options are:
#' \code{c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA")}
#'
#' \describe{
#'     \item{DCA}{Performs detrended correspondence analysis using\code{\link{decorana}}}
#'     \item{CCA}{Performs correspondence analysis, 
#'           or optionally, constrained correspondence analysis
#'           (a.k.a. canonical correspondence analysis), 
#'           via \code{\link[vegan]{cca}}}
#'     \item{RDA}{Performs redundancy analysis, or optionally 
#'           principal components analysis, via \code{\link[vegan]{rda}}}
#'     \item{CAP}{[Partial] Constrained Analysis of Principal Coordinates 
#'           or distance-based RDA, via \code{\link[vegan]{capscale}}.
#'           See \code{\link[phyloseq]{capscale.phyloseq}} for more details.
#'           In particular, a \code{\link{formula}} argument must be provided.} 
#'     \item{DPCoA}{Performs Double Principle Coordinate Analysis using a 
#'           (corrected, if necessary) phylogenetic/patristic distance
#'           between species. The calculation is performed by 
#'           \code{\link{DPCoA}}(), which ultimately uses
#'           \code{\link[ade4]{dpcoa}} after making the appropriate 
#'           accessions/corrections of the data.}
#'     \item{NMDS}{Performs Non-metric MultiDimenstional Scaling of a sample-wise 
#'           ecological distance matrix onto a user-specified number of axes, \code{k}.
#'           By default, \code{k=2}, but this can be modified as a supplementary argument.
#'           This method is ultimately carried out by \code{\link{metaMDS}} after the 
#'           appropriate accessions and distance calculations.
#'           Because \code{metaMDS} includes its own distance 
#'           calculation wrappers to \code{\link[vegan]{vegdist}}, and these provide
#'           additional functionality in the form of species scores,
#'           \code{ordinate} will pass-on the \code{distance} 
#'           argument to \code{metaMDS} if it is among the 
#'           supported \code{vegdist} methods. However, all distance methods
#'           supported by \code{\link{distance}} are supported here,
#'           including \code{"unifrac"} (the default) and \code{"DPCoA"}.}
#'     \item{MDS/PCoA}{Performs principal coordinate analysis 
#'           (also called principle coordinate decomposition, 
#'           multidimensional scaling (MDS), or classical scaling)
#'           of a distance matrix (Gower 1966), 
#'           including two correction methods for negative eigenvalues.
#'           See 
#'           \code{\link[ape]{pcoa}} for further details.
#'          }	
#'	}
#' 
#' @param distance (Optional). A character string. Default is \code{"bray"}.
#'  The name of a supported \code{\link{distance}} method; 
#'  or, alternatively, 
#'  a pre-computed \code{\link{dist}}-class object.
#'  This argument is only utilized
#'  if a distance matrix is required by the ordination method specified by the
#'  \code{method} argument (above).
#'
#'  Any supported \code{\link{distance}} methods 
#'  are supported arguments to \code{distance} here. 
#'  See \code{\link{distance}} for more details, examples.
#' 
#' @param formula (Optional). A model \code{\link{formula}}.
#'  Only relevant for certain ordination methods.
#'  The left hand side is ignored, defined by 
#'  the \code{physeq} and \code{distance} arguemnts.
#'  The right hand side gives the constraining variables,
#'  and conditioning variables can be given 
#'  within a special function \code{Condition}.
#'  See \code{\link[vegan]{cca}} or \code{\link[vegan]{capscale}}
#'  for examples/details.
#' 
#' @param ... (Optional). Additional arguments to supporting functions. For 
#'  example, the additional argument \code{weighted=TRUE} would be passed on
#'  to \code{\link{UniFrac}} if \code{"unifrac"} were chosen as the
#'  \code{distance} option and \code{"MDS"} as the ordination \code{method}
#'  option. Alternatively, if \code{"DCA"} were chosen as the 
#'  ordination \code{method} option, additional arguments would be passed on
#'  to the relevant ordination function, \code{\link{decorana}}, for example.
#' 
#' @return
#'  An ordination object. The specific class of the returned object depends upon the 
#'  ordination method, as well as the function/package that is called internally
#'  to perform it.
#'  As a general rule, any of the ordination classes
#'  returned by this function will be recognized by downstream tools in the
#'  \code{phyloseq} package, for example the ordination plotting 
#'  function, \code{\link{plot_ordination}}.
#' 
#' @seealso 
#'  \href{http://joey711.github.io/phyloseq/plot_ordination-examples}{The plot_ordination Tutorial}
#' 
#'  Related component ordination functions described within phyloseq:
#'
#'  \code{\link{DPCoA}}
#'
#'  Described/provided by other packages:
#' 
#'  \code{\link{cca}}/\code{\link{rda}}, \code{\link{decorana}}, \code{\link{metaMDS}}, 
#'  \code{\link{pcoa}}, \code{\link[vegan]{capscale}}
#'
#'  NMDS and MDS/PCoA both operate on distance matrices, typically based on some
#'  pairwise comparison of the microbiomes in an experiment/project. There are
#'  a number of common methods to use to calculate these pairwise distances, and
#'  the most convenient function (from a \code{phyloseq} point of view) for calculating 
#'  these distance matrices is the 
#'
#'  \code{\link{distance}}
#'
#'  function. It can be 
#'  thought of as a distance / dissimilarity-index companion function for
#'  \code{ordinate}, and indeed the distance options provided to \code{ordinate}
#'  are often simply passed on to \code{\link{distance}}.
#'
#'  A good quick summary of ordination is provided in the introductory vignette
#'  for vegan:
#' 
#'  \href{http://cran.r-project.org/web/packages/vegan/vignettes/intro-vegan.pdf}{vegan introductory vignette}
#'
#'  The following \code{R} task views are also useful for understanding the 
#'  available tools in \code{R}:
#'
#' \href{http://cran.r-project.org/web/views/Environmetrics.html}{Analysis of Ecological and Environmental Data}
#'
#' \href{http://cran.r-project.org/web/views/Multivariate.html}{Multivariate Statistics}
#' 
#' @importFrom vegan decorana
#' @importFrom vegan metaMDS
#' @importFrom vegan wisconsin
#' @importFrom vegan decostand
#' @importFrom ape pcoa 
#' @export
#' @examples
#' # See http://joey711.github.io/phyloseq/plot_ordination-examples
#' # for many more examples.
#' # plot_ordination(GP, ordinate(GP, "DCA"), "samples", color="SampleType")
ordinate = function(physeq, method="DCA", distance="bray", formula=NULL, ...){
  # If `physeq` is a formula, post deprecated notice, attempt to convert and dispatch
  if( inherits(physeq, "formula") ){
    .Deprecated(msg=paste0("First argument, `physeq`, as formula is deprecated.\n",
                           "There is now an explicit `formula` argument.\n",
                           "Please revise method call accordingly."))
    # Create the new formula, RHS-only
    formchar = as.character(physeq)
    # Error if only RHS. Formula-first syntax required both sides.
    if(length(formchar) < 3){
      stop("Need both sides of formula in this deprecated syntax... Revisit ordinate() documentation / examples.")
    }
    # Replace with (presumed) phyloseq object.
    physeq <- get(as.character(physeq)[2])
    # Create the new formula, RHS-only. 
    newFormula = as.formula(paste0("~", formchar[length(formchar)]))  
    # Dispatch to (hopefully) ordinate,phyloseq
    return(ordinate(physeq, method=method, distance=distance, formula=newFormula, ...))
  }
	# Define table of currently-supported methods
	method_table <- c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA")
	# List supported method names to user, if requested.
	if( inherits(physeq, "character") ){
		if( physeq=="help" ){
			cat("Available arguments to methods:\n")
			print(c(method_table))
			cat("Please be exact, partial-matching not supported.\n")
			cat("Can alternatively provide a custom distance.\n")
			cat("See:\n help(\"distance\") \n")
			return()
		} else if( physeq=="list" ){
			return(c(method_table))
		} else {
			cat("physeq needs to be a phyloseq-class object, \n")
			cat("or a character string matching \"help\" or \"list\". \n")			
		}	
	}
  # Final check that `physeq` is a phyloseq or otu_table class
  if( !inherits(physeq, "phyloseq") & !inherits(physeq, "otu_table") ){
    stop("Expected a phyloseq object or otu_table object.")
  }
	# # Start with methods that don't require 
	# #  additional distance calculation. (distance argument ignored)
	# DCA
	if( method == "DCA" ){
		return( decorana(veganifyOTU(physeq), ...) )
	}
	# CCA / RDA
	if( method %in% c("CCA", "RDA") ){
	  return(cca.phyloseq(physeq, formula, method, ...))
	}
	# CAP
	if( method == "CAP" ){
    # Call/return with do.call
	  return(capscale.phyloseq(physeq, formula, distance, ...))
	}
	# DPCoA
	if( method == "DPCoA" ){
		return( DPCoA(physeq, ...) )
	}  
	# # Now resort to methods that do require a separate distance/dist-calc
	# Define ps.dist. Check the class of distance argument is character or dist
	if( inherits(distance, "dist") ){
		ps.dist <- distance
	} else if( class(distance) == "character" ){
		# There are some special options for NMDS/metaMDS if distance-method
		# is supported by vegdist, so check first. If not, just calculate distance	
		vegdist_methods <- c("manhattan", "euclidean", "canberra", "bray", 
		"kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", 
		"mountford", "raup" , "binomial", "chao")			
		# NMDS with vegdist-method to include species		
		if(method == "NMDS" & distance %in% vegdist_methods){
			return(metaMDS(veganifyOTU(physeq), distance, ...))
		}
		# Calculate distance with handoff to distance()
		ps.dist <- distance(physeq, distance, ...)
	}
	# Vanilla MDS/PCoA
	if( method %in% c("PCoA", "MDS")){
		return(pcoa(ps.dist))
	}
	# NMDS with non-vegdist-method
	if(method == "NMDS"){
		return(metaMDS(ps.dist))
	}	
}
################################################################################
#' Calculate Double Principle Coordinate Analysis (DPCoA) 
#' using phylogenetic distance
#'
#' Function uses abundance (\code{\link{otu_table-class}}) and 
#' phylogenetic (\code{\link[ape]{phylo}}) components of a 
#' \code{\link{phyloseq-class}} experiment-level object
#' to perform a
#' Double Principle Coordinate Analysis (DPCoA), relying heavily on 
#' the underlying (and more general) function, \code{\link[ade4]{dpcoa}}.
#' The distance object ultimately provided is the square root of the 
#' cophenetic/patristic (\code{\link[ape]{cophenetic.phylo}}) distance 
#' between the species, which is always Euclidean. 
#'  
#' Although this distance is Euclidean, for numerical reasons it 
#' will sometimes look non-Euclidean, and a correction will be performed. 
#' See \code{correction} argument. 
#'
#' @param physeq (Required). A \code{\link{phyloseq-class}} object
#'  containing, at a minimum, abundance (\code{\link{otu_table-class}}) and 
#'  phylogenetic (\code{\link[ape]{phylo}}) components.
#'  As a test, the accessors \code{\link{otu_table}} and \code{\link{phy_tree}}
#'  should return an object without error.
#' 
#' @param correction (Optional). A function. The function must be
#'  able to take a non-Euclidean \code{\link{dist}}ance object,
#'  and return a new \code{dist}ance object that is Euclidean.
#'  If testing a distance object, try \code{\link[ade4]{is.euclid}}.
#' 
#'  
#'  Although the distance matrix should always be Euclidean, for numerical 
#' reasons it will sometimes appear non-Euclidean and a correction method must 
#' be applied. Two recommended correction methods are
#'  \code{\link[ade4]{cailliez}} and \code{\link[ade4]{lingoes}}.
#'  The default is \code{cailliez},
#'  but not for any particularly special reason. If the
#'  distance matrix is Euclidian, no correction will be 
#'  performed, regardless of the value of the \code{correction} argument.
#'
#' @param scannf (Optional). Logical. Default is \code{FALSE}. This
#'  is passed directly to \code{\link[ade4]{dpcoa}}, and causes a
#'  barplot of eigenvalues to be created if \code{TRUE}. This is not
#'  included in \code{...} because the default for \code{\link[ade4]{dpcoa}}
#'  is \code{TRUE}, although in many expected situations we would want
#'  to suppress creating the barplot.
#'
#' @param ... Additional arguments passed to \code{\link[ade4]{dpcoa}}.
#'
#' @return A \code{dpcoa}-class object (see \code{\link[ade4]{dpcoa}}).
#'
#' @seealso \code{\link[ade4]{dpcoa}}
#'
#' @author Julia Fukuyama \email{julia.fukuyama@@gmail.com}.
#'  Adapted for phyloseq by Paul J. McMurdie.
#'
#' @importFrom ape cophenetic.phylo
#' @importFrom ade4 cailliez
#' @importFrom ade4 dpcoa
#' @importFrom ade4 is.euclid
#' @export
#' @references
#' Pavoine, S., Dufour, A.B. and Chessel, D. (2004) 
#' From dissimilarities among species to dissimilarities among communities: 
#' a double principal coordinate analysis.
#' Journal of Theoretical Biology, 228, 523-537.
#' 
#' @examples
#' # # # # # # Esophagus
#' data(esophagus)
#' eso.dpcoa <- DPCoA(esophagus)
#' eso.dpcoa
#' plot_ordination(esophagus, eso.dpcoa, "samples")
#' plot_ordination(esophagus, eso.dpcoa, "species")
#' plot_ordination(esophagus, eso.dpcoa, "biplot")
#' #
#' #
#' # # # # # # GlobalPatterns
#' data(GlobalPatterns)
#' # subset GP to top-150 taxa (to save computation time in example)
#' keepTaxa <- names(sort(taxa_sums(GlobalPatterns), TRUE)[1:150])
#' GP       <- prune_taxa(keepTaxa, GlobalPatterns)
#' # Perform DPCoA
#' GP.dpcoa <- DPCoA(GP)
#' plot_ordination(GP, GP.dpcoa, color="SampleType")
DPCoA <- function(physeq, correction=cailliez, scannf=FALSE, ...){
	# Check that physeq is a phyloseq-class
	if(!class(physeq)=="phyloseq"){stop("physeq must be phyloseq-class")}
	
	# Remove any OTUs that are absent from all the samples.
	physeq <- prune_taxa((taxa_sums(physeq) > 0), physeq)
	
	# Access components for handing-off
	OTU  <- otu_table(physeq)
	tree <- phy_tree(physeq)
	
	# Enforce that OTU is in samples-by-species orientation
	if(taxa_are_rows(OTU) ){ OTU <- t(OTU) }
  
	# get the patristic distances between the species from the tree 
	patristicDist <- sqrt(as.dist(cophenetic.phylo(tree)))
	
	# if the patristic distances are not Euclidean, 
	# then correct them or throw meaningful error.
	if( !is.euclid(patristicDist) ){
		patristicDist <- correction(patristicDist)
		
		# Check that this is now Euclidean.
		if( !is.euclid(patristicDist) ){
			stop('Corrected distance still not Euclidean \n',
			"please provide a different correction method")
		}
	}
	
	# NOTE: the dpcoa function in ade4 requires a data.frame
	return( dpcoa(data.frame(OTU), patristicDist, scannf, ...) )
}
################################################################################
################################################################################
# vegan::cca "extension".
# formula is main input to this function. This complicates signature handling.
# A new method with a separate name is defined instead.
#
# Must transpose the phyloseq otu_table to fit the vegan::cca convention
# Whether-or-not to transpose needs to be a check, based on the 
#   "taxa_are_rows" slot value
################################################################################
#' Constrained Correspondence Analysis and Redundancy Analysis.
#' 
#' This is the internal function that simplifies getting phyloseq data
#' into the constrained ordination functions, 
#' \code{\link[vegan]{cca}} and \code{\link[vegan]{rda}}.
#' Unlike \code{\link[phyloseq]{capscale.phyloseq}}, the formula argument
#' to these methods is optional, and results in an unconstrained ordination.
#' 
#' @param physeq (Required). Phylogenetic sequencing data
#'  (\code{\link{phyloseq-class}}).
#'  The data on which you want to perform the ordination. 
#' 
#' @param formula (Optional). A \code{\link{formula}}, 
#'  specifying the contraining variable(s) format,
#'  with variable names corresponding to \code{\link{sample_data}} (RHS)
#'  from within \code{physeq}.
#'  
#' @param method (Optional). A single \code{\link{character}} string,
#' specifying \code{"RDA"} or \code{"CCA"}. Default is \code{"CCA"}.
#'
#' @param ... (Optional). Additional named arguments passed to 
#'  \code{\link[vegan]{capscale}}.
#'
#' @return same output as \code{\link[vegan]{cca}} 
#'  or \code{\link[vegan]{rda}}, respectively.
#'
#' @seealso \code{\link{plot_ordination}},
#'  \code{\link[vegan]{rda}}, \code{\link[vegan]{cca}}
#'
#' @aliases cca.phyloseq rda.phyloseq
#' @rdname cca-rda-phyloseq-methods
#' @docType methods
#'
#' @keywords internal
#' @examples #
#' # cca.phyloseq(physeq, formula, method, ...)
setGeneric("cca.phyloseq", function(physeq, formula=NULL, method="CCA", ...){
  standardGeneric("cca.phyloseq")
})
#' @importFrom vegan cca
#' @importFrom vegan rda
#' @aliases cca.phyloseq,phyloseq,formula-method
#' @rdname cca-rda-phyloseq-methods
setMethod("cca.phyloseq", signature=c("phyloseq", "formula"), 
function(physeq, formula, method="CCA", ...){
  data = data.frame(sample_data(physeq, FALSE), stringsAsFactors=FALSE)
  if( length(data) < 1 ){
    stop("`physeq` argument must include non-empty `sample_data`")
  }
	OTU = veganifyOTU(physeq)
	# Create new formula. Left-hand side is ignored.
	formchar = as.character(formula)
	newFormula = as.formula(paste0("OTU ~ ", formchar[length(formchar)]))
	# Note that ade4 also has a conflicting "cca" function.
  # You don't import ade4::cca to avoid the conflict.
  if(method=="CCA"){
    return(cca(newFormula, data=data))
  } else if(method=="RDA"){
    return(rda(newFormula, data=data))
  } else {
    warning("Unsupported `method` argument. Must be 'RDA' or 'CCA'")
    return(NULL)
  }
})
#' @importFrom vegan cca
#' @aliases cca.phyloseq,otu_table-method
#' @rdname cca-rda-phyloseq-methods
setMethod("cca.phyloseq", signature="otu_table",
          function(physeq, formula=NULL, method="CCA", ...){
  # OTU table by itself indicates an unconstrained ordination is requested.
  # Formula argument is ignored.
	if(method=="CCA"){
	  return(cca(veganifyOTU(physeq)))
	} else if(method=="RDA"){
	  return(rda(veganifyOTU(physeq)))
	} else {
	  warning("Unsupported `method` argument. Must be 'RDA' or 'CCA'")
    return(NULL)
	}
})
#' @importFrom vegan cca
#' @aliases cca.phyloseq,phyloseq,NULL-method
#' @rdname cca-rda-phyloseq-methods
setMethod("cca.phyloseq", signature=c("phyloseq", "NULL"), 
function(physeq, formula, method="CCA", ...){
  # Absence of a formula (NULL) indicates unconstrained ordination.
  # Access otu_table, and dispatch.
  return(cca.phyloseq(otu_table(physeq), NULL, method, ...))
})
################################################################################
#' Estimate the gap statistic on an ordination result
#' 
#' This is a wrapper for the \code{\link[cluster]{clusGap}} function,
#' expecting an ordination result as the main data argument.
#' 
#' @param ord (Required). An ordination object. The precise class can vary.
#'  Any ordination classes supported internally by the phyloseq package 
#'  should work, ultimately by passing to the \code{\link[vegan]{scores}} function
#'  or its internal extensions in phyloseq.
#' @param axes (Optional). The ordination axes that you want to include.
#' @param type (Optional). One of \code{"sites"} 
#'  (the vegan package label for samples) or 
#'  \code{"species"} (the vegan package label for OTUs/taxa).
#'  Default is \code{"sites"}.
#' @param FUNcluster (Optional). This is passed to \code{\link[cluster]{clusGap}}.
#'  The documentation is copied here for convenience: 
#'  a function which accepts as first argument a (data) matrix like \code{x}, 
#'  second argument, say (the number of desired clusters) \code{k}, where \code{k >= 2},
#'  and returns a list with a component named (or shortened to) cluster
#'  which is a vector of length \code{n = nrow(x)} of integers in \code{1:k}
#'  determining the clustering or grouping of the \code{n} observations.
#'  The default value is the following function, which wraps
#'  partitioning around medoids, \code{\link[cluster]{pam}}:
#'  
#'  \code{function(x, k){list(cluster = pam(x, k, cluster.only=TRUE))}}
#'  
#'  Any function that has these input/output properties (performing a clustering)
#'  will suffice. The more appropriate the clustering method, the better chance
#'  your gap statistic results will be useful.
#' @param K.max	(Optional). A single positive integer value.
#'  It indicates the maximum number of clusters that will be considered.
#'  Value must be at least two.
#'  This is passed to \code{\link[cluster]{clusGap}}.
#' @param ... (Optional). Additional named parameters
#'  passed on to \code{\link[cluster]{clusGap}}.
#'  For example, the \code{method} argument provides for extensive options
#'  regarding the method by which the ``optimal'' number of clusters
#'  is computed from the gap statistics (and their standard deviations).
#'  See the \code{\link[cluster]{clusGap}} documentation for more details.
#' 
#' @return
#' An object of S3 class \code{"clusGap"}, basically a list with components.
#' See the \code{\link[cluster]{clusGap}} documentation for more details.
#' 
#' @importFrom vegan scores
#' @importFrom cluster clusGap
#' @importFrom cluster pam
#' @export
#' @examples
#' data("soilrep")
#' sord  = ordinate(soilrep, "PCoA", "bray")
#' # Evaluate axes with scree plot
#' plot_scree(sord)
#' # Gap Statistic
#' gs = gapstat_ord(sord, axes=1:3, verbose=FALSE)
#' # plot_ordination(soilrep, sord,  color="Treatment")
#' plot_clusgap(gs)
#' print(gs, method="Tibs2001SEmax")
gapstat_ord = function(ord, axes=c(1:2), type="sites", 
	FUNcluster=function(x, k){list(cluster = pam(x, k, cluster.only=TRUE))},
	K.max=8, ...){
	#
	# Use the scores function to get the ordination coordinates
	x = scores(ord, display=type)
	# If axes not explicitly defined (NULL), then use all of them
	if(is.null(axes)){
		axes = 1:ncol(x)
	}
	# Finally, perform, and return, the gap statistic calculation using
	# cluster::clusGap
	return(clusGap(x[, axes], FUNcluster, K.max, ...))
}
################################################################################
# Define an internal function for accessing and orienting the OTU table
# in a fashion suitable for vegan functions
# @keywords internal
veganifyOTU <- function(physeq){
  if(taxa_are_rows(physeq)){physeq <- t(physeq)}
  return(as(otu_table(physeq), "matrix"))
}
################################################################################
#' Constrained Analysis of Principal Coordinates, \code{\link[vegan]{capscale}}.
#'
#' See \code{\link[vegan]{capscale}} for details. A formula is main input.
#' 
#' @param physeq (Required). Phylogenetic sequencing data
#'  (\code{\link{phyloseq-class}}).
#'  The data on which you want to perform the ordination. 
#' 
#' @param formula (Required). A \code{\link{formula}}, specifying the input.
#'  No need to directly access components. \code{capscale.phyloseq} understands 
#'  where to find the abundance table (LHS) and \code{\link{sample_data}} (RHS)
#'  from within the phyloseq object. 
#'
#' @param distance (Required). A \code{\link{character}} string, specifying
#'  the name of the dissimilarity (or distance) method supported by
#'  the phyloseq \code{\link[phyloseq]{distance}} function.
#'  Alternatively, a pre-computed \code{\link{dist}}-object can be provided here,
#'  in which case it supersedes any use of the \code{\link{otu_table}}
#'  in your phyloseq object.
#'  
#'  Note that \code{\link[vegan]{capscale}} 
#'  with Euclidean distances will be identical to \code{\link[vegan]{rda}}
#'  in eigenvalues and in site, species, and biplot scores
#'  (except for possible sign reversal). However, it makes no sense to use 
#'  \code{\link[vegan]{capscale}} with Euclidean distances, 
#'  since direct use of \code{\link[vegan]{rda}} is much more efficient
#'  (and supported in the \code{\link{ordinate}} function with \code{method=="RDA"})
#'  Even with non-Euclidean dissimilarities,
#'  the rest of the analysis will be metric and linear.
#'
#' @param ... (Optional). Additional named arguments passed to 
#'  \code{\link[vegan]{capscale}}.
#'
#' @return Ordination object defined by \code{\link[vegan]{capscale}}.
#'
#' @seealso
#'  \code{\link{plot_ordination}}
#' 
#'  \code{\link[vegan]{rda}}
#'  
#'  \code{\link[vegan]{capscale}}
#'
#' @aliases capscale.phyloseq
#' @rdname capscale-phyloseq-methods
#' @docType methods
#' @importFrom vegan capscale
#' @keywords internal
#' @examples 
#' # See other examples at
#' # http://joey711.github.io/phyloseq/plot_ordination-examples
#' data(GlobalPatterns)
#' GP = prune_taxa(names(sort(taxa_sums(GlobalPatterns), TRUE)[1:50]), GlobalPatterns)
#' ordcap = ordinate(GP, "CAP", "bray", ~SampleType)
#' plot_ordination(GP, ordcap, "samples", color="SampleType")
setGeneric("capscale.phyloseq", function(physeq, formula, distance, ...){
  data = data.frame(sample_data(physeq, FALSE), stringsAsFactors=FALSE)
  if( length(data) < 1 ){
    stop("`physeq` argument must include non-empty `sample_data`")
  }
  standardGeneric("capscale.phyloseq")
})
#' @importFrom vegan capscale
#' @aliases capscale.phyloseq,phyloseq,formula,dist-method
#' @rdname capscale-phyloseq-methods
setMethod("capscale.phyloseq", c("phyloseq", "formula", "dist"),
function(physeq, formula, distance, ...){
  data = data.frame(sample_data(physeq), stringsAsFactors=FALSE)
  # Convert formula to character vector, compute on language.
  formchar = as.character(formula)  
  newFormula = as.formula(paste0("distance ~ ", formchar[length(formchar)]))
  return(capscale(formula=newFormula, data=data, ...))
})
#' @importFrom vegan capscale
#' @aliases capscale.phyloseq,phyloseq,formula,character-method
#' @rdname capscale-phyloseq-methods
setMethod("capscale.phyloseq", c("phyloseq", "formula", "character"),
function(physeq, formula, distance, ...){
  data = data.frame(sample_data(physeq), stringsAsFactors=FALSE)
  # The goal here is to process the distance identifier string
  # and dispatch accordingly.
  if( length(distance) != 1 ){
    warning("`distance` was unexpected length. \n",
            " `distance` argument should be a single character string",
            " or dist matrix. \n",
            "Attempting to use first element only.")
  }
  distance <- distance[1]
  if(!distance %in% unlist(distanceMethodList)){
    # distance must be among the supported distance options
    # (which is a superset of vegdist).
    stop("The distance method you specified is not supported by phyloseq")
  }
  # Convert formula to character vector, compute on language.
  formchar = as.character(formula)
  if(distance %in% distanceMethodList$vegdist){
    # If it is among the vegdist distances, pass it along to vegan::capscale
    OTU = veganifyOTU(physeq)
    newFormula = as.formula(paste0("OTU ~ ", formchar[length(formchar)]))
    return(capscale(formula=newFormula, data=data, distance=distance, ...))
  } else {
    # Else calculate the distance matrix here, and dispatch.
    distance <- distance(physeq=physeq, method=distance, type="samples")
    return(capscale.phyloseq(physeq, formula, distance, ...))
  }
})
################################################################################
