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
#' @seealso \code{\link{plot_ordination_phyloseq}}, \code{\link{calcplot}},
#'  \code{\link[vegan]{rda}}, \code{\link[vegan]{cca}}
#'
#' @aliases cca.phyloseq rda.phyloseq
#' @rdname cca-rda-phyloseq-methods
#' @docType methods
#'
#' @export
#' @import vegan
#' @examples #
#' # data(ex1)
#' # # For RDA, use thresholded-rank
#' # ex4  <- transformsamplecounts(ex1, threshrankfun(500))
#' # # RDA
#' # modr <- rda.phyloseq(ex4 ~ Diet + Gender)
#' # # CCA
#' # modc <- cca.phyloseq(ex1 ~ Diet + Gender)
#' # plot_ordination_phyloseq(modr, ex1)
#' # plot_ordination_phyloseq(modc, ex1)
#' # # Perform unconstrained ordination
#' # mod1 <- cca.phyloseq(ex1)
#' # # unconstrained plot using vegan plotting
#' # vegan:::plot.cca(mod1)
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