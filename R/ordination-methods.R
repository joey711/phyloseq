################################################################################
#' Calculate Double Principle Coordinate Analysis (DPCoA) 
#' using phylogenetic distance
#'
#' Function uses abundance (\code{\link{otuTable-class}}) and 
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
#'  containing, at a minimum, abundance (\code{\link{otuTable-class}}) and 
#'  phylogenetic (\code{\link[ape]{phylo}}) components.
#'  As a test, the accessors \code{\link{otuTable}} and \code{\link{tre}}
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
#' # keepTaxa <- names(sort(speciesSums(GlobalPatterns), TRUE)[1:150])
#' # GP       <- prune_species(keepTaxa, GlobalPatterns)
#' # # Perform DPCoA
#' # GP.dpcoa <- DPCoA(GP)
#' # plot_ordination(GP, GP.dpcoa, color="SampleType")
DPCoA <- function(physeq, correction=cailliez, scannf=FALSE, ...){
	# Check that physeq is a phyloseq-class
	if(!class(physeq)=="phyloseq"){stop("physeq must be phyloseq-class")}
	
	# Remove any OTUs that are absent from all the samples.
	physeq <- prune_species((speciesSums(physeq) > 0), physeq)
	
	# Access components for handing-off
	OTU  <- otuTable(physeq)
	tree <- tre(physeq)
	
	# Enforce that OTU is in species-by-samples orientation
	if( !speciesAreRows(OTU) ){ OTU <- t(OTU) }
  
	# get the patristic distances between the species from the tree 
	patristicDist <- as.dist(ape::cophenetic.phylo(tree))
	
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
# Must transpose the phyloseq otuTable to fit the vegan::cca convention
# Whether-or-not to transpose needs to be a check, based on the 
#   "SpeciesAreRows" slot value
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
#'  \code{\link{otuTable-class}} or \code{\link{phyloseq-class}} (without
#'  the \code{~} signifying a formula), in which case an unconstrained ordination
#'  is performed. 
#'
#' @param ... (Optional). E.g. \code{data=DF}, where \code{DF} is a \code{data.frame}
#'  containing information equivalent to
#'  a \code{sampleData} object / component. Only necessary if complex object
#'  does not already contain \code{sampleData} or you are keeping the data 
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
#' @export
#' @import vegan
#' @examples #
#' # data(GlobalPatterns)
#' # # For RDA, use thresholded-rank
#' # ex4  <- transformsamplecounts(GlobalPatterns, threshrankfun(500))
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
	OTU    <- otuTable( physeq )
	if( speciesAreRows(OTU) ){
		OTU <- t(as(OTU, "matrix"))
	} else {
		OTU <- as(OTU, "matrix")
	}
	# Create the new formula
	newFormula = as.formula(paste("OTU", as.character(X)[3], sep=" ~ "))
	# If an alternative table is not provided, assume it is from the sampleData slot
	if( is.null(data) ){
		data <- data.frame(sampleData(physeq))
	}
	# Good idea to qualify, as ade4 also has a conflicting "cca"
	# and might be a dependency in the future.	
	vegan::cca(newFormula, data=data)	
})
################################################################################
#' @aliases cca.phyloseq,otuTable-method
#' @rdname cca-rda-phyloseq-methods
setMethod("cca.phyloseq", "otuTable", function(X){
	if( speciesAreRows(X) ){
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
	cca.phyloseq(otuTable(X))
})
################################################################################
#' @usage rda.phyloseq(X, ...)
#' @export
#' @import vegan
#' @rdname cca-rda-phyloseq-methods
#' @aliases cca.phyloseq rda.phyloseq
setGeneric("rda.phyloseq", function(X, ...) standardGeneric("rda.phyloseq"))
#' @aliases rda.phyloseq,formula-method
#' @rdname cca-rda-phyloseq-methods
setMethod("rda.phyloseq", "formula", function(X, data=NULL){
	#require(vegan)
	physeq <- get( as.character(X)[2] )
	OTU    <- otuTable( physeq )
	if( speciesAreRows(OTU) ){
		OTU <- as(t(OTU), "matrix")
	} else {
		OTU <- as(OTU, "matrix")
	}
	# Create the new formula
	newFormula = as.formula(paste("OTU", as.character(X)[3], sep=" ~ "))
	# If an alternative table is not provided, assume it is from the sampleData slot
	if( is.null(data) ){
		data <- data.frame(sampleData(physeq))
	}
	rda(newFormula, data=data)	
})
################################################################################
#' @aliases rda.phyloseq,otuTable-method
#' @rdname cca-rda-phyloseq-methods
setMethod("rda.phyloseq", "otuTable", function(X){
	if( speciesAreRows(X) ){
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
	rda.phyloseq(otuTable(X))
})
################################################################################