####################################################################################
#' Wrapper of mt.maxT and mt.minP for the phyloseq package.
#'
#' @param X An \code{otuTable} or an object that contains an \code{otuTable}.
#'  In this multiple testing framework, different taxa correspond to variables
#'  (hypotheses), and samples to observations.
#'
#' @param classlabel Either a single character string indicating the column
#'  in the \code{sampleMap} of X, or, alternatively, an integer vector
#'  of length equal to the number of samples of X. The latter option is
#'  equivalent to \code{classlabel} of \code{mt.maxT} in \code{multtest}. A 
#'  third option is for classlabel to be a 2-level factor, with length
#'  equal to the number of samples of X. 
#'
#' @param ... Additional arguments, forwarded to \code{mt.maxT} 
#' 
#' @return A dataframe with components specified in the documentation for
#'  \code{mt.maxT} or \code{mt.minP} respectively.
#'
#' @seealso mt.maxT mt.minP
#'
#' @import multtest
#' @export
#' @examples #
#' ## data(ex1)
#' ## mt(ex1, "Diet")
#' ## mt(phyloseqTax(ex1), "Diet", minPmaxT="minsldzxc")
setGeneric("mt", function(X, classlabel, minPmaxT="minP", ...)
	standardGeneric("mt")
)
# otuTree doesn't have a sampleMap. just grab otuTable hope 4 best.
setMethod("mt", "otuTree", function(X, classlabel, minPmaxT="minP", ...){
	mt(otuTable(X), classlabel, minPmaxT, ...)
})
# phyloseqTree does have a sampleMap. dispatch it to phyloseq version.
setMethod("mt", "phyloseqTree", function(X, classlabel, minPmaxT="minP", ...){
	mt(phyloseq(X), classlabel, minPmaxT, ...)
})
# phyloseqTree does have a sampleMap. dispatch it to phyloseq version.
setMethod("mt", "phyloseqTax", function(X, classlabel, minPmaxT="minP", ...){
	mt(phyloseq(X), classlabel, minPmaxT, ...)
})
setMethod("mt", c("phyloseq", "factor"),
	function(X, classlabel, minPmaxT="minP", ...){
	mt(otuTable(X), classlabel, minPmaxT, ...)
})
setMethod("mt", c("phyloseq", "integer"),
	function(X, classlabel, minPmaxT="minP", ...){
	mt(otuTable(X), classlabel, minPmaxT, ...)
})
setMethod("mt", c("phyloseq", "character"),
	function(X, classlabel, minPmaxT="minP", ...){
	rawFactor <- as(sampleMap(X), "data.frame")[, classlabel]
	mt(otuTable(X), factor(rawFactor), minPmaxT, ...)
})
setMethod("mt", c("otuTable", "factor"),
	function(X, classlabel, minPmaxT="minP", ...){
	# integerize classlabel to logical suitable for mtmaxT
	classlabel <- (0:1)[classlabel]
	# Use dispatch of mtmaxT with classlabel now a suitable classlabel
	mt(X, classlabel, minPmaxT, ...)
})	
setMethod("mt", "otuTable", function(X, classlabel, minPmaxT="minP", ...){
	# require(multtest)
	if( !speciesAreRows(X) ){
		X <- t(X)
	}
	# Now call mt.maxT or mt.minP
	if( minPmaxT == "minP" ){
		mt.minP(as(X, "matrix"), classlabel, minPmaxT, ...)
	} else if( minPmaxT == "maxT" ){
		mt.maxT(as(X, "matrix"), classlabel, minPmaxT, ...)	
	} else {
		print("Nothing calculated. minPmaxT argument must be either minP or maxT.")
	}
})
####################################################################################
# #' Wrapper of mt.minP for the phyloseq package.
# #'
# #' @param X An \code{otuTable} or an object that contains an \code{otuTable}.
# #'  In this multiple testing framework, different taxa correspond to variables
# #'  (hypotheses), and samples to observations.
# #'
# #' @param classlabel Either a single character string indicating the column
# #'  in the \code{sampleMap} of X, or, alternatively, an integer vector
# #'  of length equal to the number of samples of X. The latter option is
# #'  equivalent to \code{classlabel} of \code{mt.minP} in \code{multtest}. A 
# #'  third option is for classlabel to be a 2-level factor, with length
# #'  equal to the number of samples of X. 
# #'
# #' @param ... Additional arguments, forwarded to \code{mt.minP} 
# #' 
# #' @seealso mtmaxT mt.maxT
# #'
# #' @import multtest
# #' @export
# #' @examples #
# #' #data(ex1)
# #' #mtminP(ex1, "Diet")
# setGeneric("mtminP", function(X, classlabel, ...) standardGeneric("mtminP"))
# setMethod("mtminP", "otuTree", function(X, classlabel, ...){
	# mtminP(otuTable(X), classlabel, ...)
# })
# setMethod("mtminP", "phyloseq", function(X, classlabel, ...){
	# mtminP(otuTable(X), classlabel, ...)	
# })
# setMethod("mtminP", c("otuTable", "character"), function(X, classlabel, ...){
	# rawFactor <- as(sampleMap(X), "matrix")[, classlabel]
	# # Use dispatch of classlabel as factor
	# mtminP(X, factor(rawFactor), ...)
# })
# setMethod("mtminP", c("otuTable", "factor"), function(X, classlabel, ...){
	# # integerize classlabel to logical suitable for mtminP
	# classlabel <- (0:1)[classlabel]
	# # Use dispatch of mtminP with classlabel now a suitable classlabel
	# mtminP(X, classlabel, ...)
# })	
# setMethod("mtminP", "otuTable", function(X, classlabel, ...){
	# # require(multtest)
	# if( !speciesAreRows(X) ){
		# otuTable(X) <- t(otuTable(X))
	# }
	# # Now call mt.minP
	# mt.minP(as(X, "matrix"), classlabel, ...)
# })
# ####################################################################################
