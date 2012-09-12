################################################################################
#' Perform an ordination on phyloseq data
#'
#' This function wraps several commonly-used ordination methods. The type of 
#' ordination depends upon the argument to \code{method}. Try 
#' \code{ordinate("help")} or \code{ordinate("list")} for the currently
#' supported method options.
#'
#' @usage ordinate(physeq, method="DCA", distance="unifrac", ...)
#'
#' @param physeq (Required). Phylogenetic sequencing data
#'  (\code{\link{phyloseq-class}}). The data on which you want to perform the
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
#' \code{c("DCA", "CCA", "RDA", "DPCoA", "NMDS", "MDS", "PCoA")}
#'
#' \describe{
#'     \item{DCA}{Performs detrended correspondence analysis using \code{\link{decorana}}}
#'     \item{CCA}{Performs correspondence analysis, or optionally, constrained correspondence
#'           analysis (a.k.a. canonical correspondence analysis), via \code{\link[vegan]{cca}}
#'          }
#'     \item{RDA}{Performs redundancy analysis, or optionally principal components analysis,
#'           via \code{\link[vegan]{rda}}
#'          }
#'     \item{DPCoA}{Performs Double Principle Coordinate Analysis using a 
#'           (corrected, if necessary) phylogenetic/patristic distance between species. The
#'           calculation is performed by \code{\link{DPCoA}}(), which ultimately uses
#'           \code{\link[ade4]{dpcoa}} after making the appropriate 
#'           accessions/corrections of the data.
#'          }
#'     \item{NMDS}{Performs Non-metric MultiDimenstional Scaling of a sample-wise 
#'           ecological distance matrix onto a user-specified number of axes, \code{k}.
#'           By default, \code{k=2}, but this can be modified as a supplementary argument.
#'           This method is ultimately carried out by \code{\link{metaMDS}} after the 
#'           appropriate accessions and distance calculations. Because \code{metaMDS} includes
#'           its own distance calculation wrappers to \code{\link[vegan]{vegdist}}, and these provide
#'           additional functionality in the form of species scores, \code{ordinate} will
#'           pass-on the \code{distance} argument to \code{metaMDS} if it is among the 
#'           supported \code{vegdist} methods. However, all distance methods supported by
#'           \code{\link{distance}} are supported here, including \code{"unifrac"}
#'           (the default) and \code{"DPCoA"}.
#'          }
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
#' @param distance (Optional). A character string matching a 
#'  \code{\link{distance}} method; or, alternatively, 
#'  a pre-computed \code{\link{dist}}-class object.
#'  This argument is only utilized
#'  if a distance matrix is required by the ordination method specified by the
#'  \code{method} argument (above).
#'
#'  Any supported \code{\link{distance}} methods 
#'  are supported arguments to \code{distance} here. 
#'  Try \code{distance("list")} for a explicitly supported distance method
#'  abbreviations. User-specified custom distance equations should also work,
#'  e.g. \code{"(A+B-2*J)/(A+B)"}. 
#'  See \code{\link{distance}} for more details, examples.
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
#'  Related component ordination functions described within phyloseq:
#'
#'  \code{\link{DPCoA}}
#'
#'  Described/provided by other packages:
#' 
#'  \code{\link{cca}}/\code{\link{rda}}, \code{\link{decorana}}, \code{\link{metaMDS}}, 
#'  \code{\link{pcoa}}
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
#'  simply passed on to \code{\link{distance}}.
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
#' @import vegan
#' @import ape
#' @export
#' @examples
#' # # Take a subset of the GP dataset for quicker computation of examples
#' # data(GlobalPatterns)
#' # # Keep top 200 species
#' # topsp <- names(sort(taxa_sums(GlobalPatterns), TRUE)[1:200])
#' # GP    <- prune_taxa(topsp, GlobalPatterns)
#' # # Subset further to top 5 phyla
#' # top5ph <- sort(tapply(taxa_sums(GP), tax_table(GP)[, "Phylum"], sum), decreasing=TRUE)[1:5]
#' # GP     <- subset_taxa(GP, Phylum %in% names(top5ph)) 
#' # # 
#' # # Examples performing ordination with NMDS. Default distance is unweighted UniFrac
#' # GP.NMDS.UF.ord   <- ordinate(GP, "NMDS")
#' # GP.NMDS.wUF.ord  <- ordinate(GP, "NMDS", "unifrac", weighted=TRUE)
#' # GP.NMDS.Bray.ord <- ordinate(GP, "NMDS", "bray")
#' # # 
#' # # # An example plot with default, or manually-defined shapes
#' # (p <- plot_ordination(GP, GP.NMDS.Bray.ord, "biplot", color="SampleType", shape="Phylum"))
#' # # define manual shape scale:
#' # man.shapes <- 21:25
#' # names(man.shapes) <- c(get_taxa_unique(GP, "Phylum"))
#' # man.shapes <- c(samples=19, man.shapes)
#' # p + scale_shape_manual(value=man.shapes)
#' # # 
#' # # An example of constrained ordination
#' # GP.cca <- ordinate(GP~SampleType, "CCA")
#' # # 
#' # # Run-through "quick" plot examples of the other ordination options currently supported
#' # # Only showing "samples" in these examples, but "species" options supported for some methods
#' # plot_ordination(GP, ordinate(GP, "DCA"), "samples", color="SampleType")
#' # plot_ordination(GP, ordinate(GP, "CCA"), "samples", color="SampleType")
#' # plot_ordination(GP, ordinate(GP~SampleType, "CCA"), "samples", color="SampleType")
#' # plot_ordination(GP, ordinate(GP, "RDA"), "samples", color="SampleType")
#' # plot_ordination(GP, ordinate(GP~SampleType, "RDA"), "samples", color="SampleType")
#' # plot_ordination(GP, ordinate(GP, "DPCoA"), "samples", color="SampleType")
#' # plot_ordination(GP, ordinate(GP, "MDS"), "samples", color="SampleType")
#' # plot_ordination(GP, ordinate(GP, "PCoA"), "samples", color="SampleType")
#' # plot_ordination(GP, ordinate(GP, "NMDS"), "samples", color="SampleType")
#' # plot_ordination(GP, ordinate(GP, "NMDS", "w"), "samples", color="SampleType")
ordinate <- function(physeq, method="DCA", distance="unifrac", ...){
	# The table of currently-supported methods
	method_table <- c("DCA", "CCA", "RDA", "DPCoA", "NMDS", "MDS", "PCoA")
	# List supported method names to user, if requested.
	if(class(physeq) == "character"){
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

	# Define an internal function for accessing and orienting the OTU table
	# in a fashion suitable for vegan/picante functions
	veganify <- function(physeq){
		OTU <- otu_table(physeq)
		if( taxa_are_rows(OTU) ){ OTU <- t(OTU) }
		return( as(OTU, "matrix") )
	}

	# # Start with methods that don't require 
	# #  additional distance calculation. (distance argument ignored)
	# DCA
	if( method == "DCA" ){
		return( decorana(veganify(physeq), ...) )
	}
	# CCA
	if( method == "CCA" ){
		return( cca.phyloseq(physeq, ...) )
	}
	# RDA
	if( method == "RDA" ){
		return( rda.phyloseq(physeq, ...) )
	}
	# DPCoA
	if( method == "DPCoA" ){
		return( DPCoA(physeq, ...) )
	}	

	# # Now resort to methods that do require a separate distance/dist-calc
	# Define ps.dist. Check the class of distance argument is character or dist
	if( class(distance) == "dist" ){
		ps.dist <- distance
	} else if( class(distance) == "character" ){
		# There are some special options for NMDS/metaMDS if distance-method
		# is supported by vegdist, so check first. If not, just calculate distance	
		vegdist_methods <- c("manhattan", "euclidean", "canberra", "bray", 
		"kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", 
		"mountford", "raup" , "binomial", "chao")			
		# NMDS with vegdist-method to include species		
		if(method == "NMDS" & distance %in% vegdist_methods){
			return(metaMDS(veganify(physeq), distance, ...))
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
#' The distance object ultimately provided as the cophenetic/patristic
#' (\code{\link[ape]{cophenetic.phylo}}) distance between the species. 
#' 
#' In most real-life, real-data applications, the phylogenetic tree 
#' will not provide a Euclidean distance matrix, and so a correction
#' will be performed, if needed. See \code{correction} argument. 
#'
#' @usage DPCoA(physeq, correction=cailliez, scannf=FALSE, ...)
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
#'  In most real-life, real-data applications, the phylogenetic tree 
#'  will not provide a Euclidean distance matrix, and so a correction
#'  will be needed. 
#'  Two recommended correction methods are
#'  \code{\link[ade4]{cailliez}} and \code{\link[ade4]{lingoes}}.
#'  The default is \code{cailliez},
#'  but not for any particularly special reason. If the patristic 
#'  distance matrix turns out to be Euclidian, no correction will be 
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
#' @import ape
#' @export
#' @references
#' Pavoine, S., Dufour, A.B. and Chessel, D. (2004) 
#' From dissimilarities among species to dissimilarities among communities: 
#' a double principal coordinate analysis.
#' Journal of Theoretical Biology, 228, 523-537.
#' 
#' @examples
#' # # # # # # # Esophagus
#' # data(esophagus)
#' # eso.dpcoa <- DPCoA(esophagus)
#' # plot_ordination(esophagus, eso.dpcoa, "samples")
#' # plot_ordination(esophagus, eso.dpcoa, "species")
#' # plot_ordination(esophagus, eso.dpcoa, "biplot")
#' # #
#' # #
#' # # # # # # # GlobalPatterns
#' # data(GlobalPatterns)
#' # # subset GP to top-150 taxa (to save computation time in example)
#' # keepTaxa <- names(sort(taxa_sums(GlobalPatterns), TRUE)[1:150])
#' # GP       <- prune_taxa(keepTaxa, GlobalPatterns)
#' # # Perform DPCoA
#' # GP.dpcoa <- DPCoA(GP)
#' # plot_ordination(GP, GP.dpcoa, color="SampleType")
DPCoA <- function(physeq, correction=cailliez, scannf=FALSE, ...){
	# Check that physeq is a phyloseq-class
	if(!class(physeq)=="phyloseq"){stop("physeq must be phyloseq-class")}
	
	# Remove any OTUs that are absent from all the samples.
	physeq <- prune_taxa((taxa_sums(physeq) > 0), physeq)
	
	# Access components for handing-off
	OTU  <- otu_table(physeq)
	tree <- phy_tree(physeq)
	
	# Enforce that OTU is in species-by-samples orientation
	if( !taxa_are_rows(OTU) ){ OTU <- t(OTU) }
  
	# get the patristic distances between the species from the tree 
	patristicDist <- as.dist(cophenetic.phylo(tree))
	
	# if the patristic distances are not Euclidean, 
	# then correct them or throw meaningful error.
	if( !ade4::is.euclid(patristicDist) ){
		patristicDist <- correction(patristicDist)
		
		# Check that this is now Euclidean.
		if( !ade4::is.euclid(patristicDist) ){
			stop('Corrected distance still not Euclidean \n',
			"please provide a different correction method")
		}
	}
	
	# NOTE: the dpcoa function in ade4 requires a data.frame
	return( ade4::dpcoa(data.frame(OTU), patristicDist, scannf, ...) )
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
#' Wrapper for \code{\link[vegan]{cca}} and \code{\link[vegan]{rda}}.
#'
#' A formula is main input to \code{\link[vegan]{cca}}. This complicates dispatch based
#' on object signature. A new method with a separate name is defined instead.
#'
#' @usage cca.phyloseq(X, ...)
#' 
#' @param X (Required). A \code{\link{formula}}, specifying the input.
#'  No need to directly access components.
#'  \code{cca.phyloseq} understands where to find the abundance table
#'  and sample data. Alternatively, \code{X} can be an 
#'  \code{\link{otu_table-class}} or \code{\link{phyloseq-class}} (without
#'  the \code{~} signifying a formula), in which case an unconstrained ordination
#'  is performed. 
#'
#' @param ... (Optional). E.g. \code{data=DF}, where \code{DF} is a \code{data.frame}
#'  containing information equivalent to
#'  a \code{sample_data} object / component. Only necessary if complex object
#'  does not already contain \code{sample_data} or you are keeping the data 
#'  separate for some reason.
#'
#' @return same output as \code{\link[vegan]{cca}} or \code{\link[vegan]{rda}}, respectively.
#'
#' @seealso \code{\link{plot_ordination}},
#'  \code{\link[vegan]{rda}}, \code{\link[vegan]{cca}}
#'
#' @aliases cca.phyloseq rda.phyloseq
#' @rdname cca-rda-phyloseq-methods
#' @docType methods
#'
#' @keywords internal
#' @import vegan
#' @examples #
#' # data(GlobalPatterns)
#' # # For RDA, use thresholded-rank
#' # ex4  <- transform_sample_counts(GlobalPatterns, threshrankfun(500))
#' # # RDA
#' # modr <- rda.phyloseq(ex4 ~ SampleType)
#' # # CCA
#' # modc <- cca.phyloseq(GlobalPatterns ~ SampleType)
#' # plot_ordination(GlobalPatterns, modr, "biplot")
#' # plot_ordination(GlobalPatterns, modc, "biplot")
#' # # Perform unconstrained ordination
#' # mod1 <- cca.phyloseq(GlobalPatterns)
setGeneric("cca.phyloseq", function(X, ...) standardGeneric("cca.phyloseq"))
################################################################################
#' @aliases cca.phyloseq,formula-method
#' @rdname cca-rda-phyloseq-methods
setMethod("cca.phyloseq", "formula", function(X, data=NULL){
	physeq <- get( as.character(X)[2] )
	OTU    <- otu_table( physeq )
	if( taxa_are_rows(OTU) ){
		OTU <- t(as(OTU, "matrix"))
	} else {
		OTU <- as(OTU, "matrix")
	}
	# Create the new formula
	newFormula = as.formula(paste("OTU", as.character(X)[3], sep=" ~ "))
	# If an alternative table is not provided, assume it is from the sample_data slot
	if( is.null(data) ){
		data <- data.frame(sample_data(physeq))
	}
	# Good idea to qualify, as ade4 also has a conflicting "cca"
	# and might be a dependency in the future.	
	vegan::cca(newFormula, data=data)	
})
################################################################################
#' @aliases cca.phyloseq,otu_table-method
#' @rdname cca-rda-phyloseq-methods
setMethod("cca.phyloseq", "otu_table", function(X){
	if( taxa_are_rows(X) ){
		X <- t(as(X, "matrix"))
	} else {
		X <- as(X, "matrix")
	}
	# Good idea to qualify, as ade4 also has a conflicting "cca"
	# and might be a dependency in the future.
	vegan::cca(X)	
})
################################################################################
#' @aliases cca.phyloseq,phyloseq-method
#' @rdname cca-rda-phyloseq-methods
setMethod("cca.phyloseq", "phyloseq", function(X){
	cca.phyloseq(otu_table(X))
})
################################################################################
#' @keywords internal
#' @usage rda.phyloseq(X, ...)
#' @import vegan
#' @rdname cca-rda-phyloseq-methods
#' @aliases cca.phyloseq rda.phyloseq
setGeneric("rda.phyloseq", function(X, ...) standardGeneric("rda.phyloseq"))
#' @aliases rda.phyloseq,formula-method
#' @rdname cca-rda-phyloseq-methods
setMethod("rda.phyloseq", "formula", function(X, data=NULL){
	#require(vegan)
	physeq <- get( as.character(X)[2] )
	OTU    <- otu_table( physeq )
	if( taxa_are_rows(OTU) ){
		OTU <- as(t(OTU), "matrix")
	} else {
		OTU <- as(OTU, "matrix")
	}
	# Create the new formula
	newFormula = as.formula(paste("OTU", as.character(X)[3], sep=" ~ "))
	# If an alternative table is not provided, assume it is from the sample_data slot
	if( is.null(data) ){
		data <- data.frame(sample_data(physeq))
	}
	rda(newFormula, data=data)	
})
################################################################################
#' @aliases rda.phyloseq,otu_table-method
#' @rdname cca-rda-phyloseq-methods
setMethod("rda.phyloseq", "otu_table", function(X){
	if( taxa_are_rows(X) ){
		X <- t(as(X, "matrix"))
	} else {
		X <- as(X, "matrix")
	}
	rda(X)	
})
################################################################################
#' @aliases rda.phyloseq,phyloseq-method
#' @rdname cca-rda-phyloseq-methods
setMethod("rda.phyloseq", "phyloseq", function(X){
	rda.phyloseq(otu_table(X))
})
################################################################################
