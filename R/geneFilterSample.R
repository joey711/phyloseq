############################################################
#' Filter OTUs sample-wise.
#' 
#' A general OTU trim function for selecting OTUs that satisfy
#'  some criteria within the distribution of each sample, and then
#'  also an additional criteria for number of samples.
#' By contrast, \code{genefilter}, of the genefilter package in Bioconductor,
#' works only on the rows of a matrix.
#' Here we want a genefilter function that considers sample-wide (column-wide)
#' criteria and determines which rows are acceptable
#' within each column. The number of acceptable samples (columns) is then used
#' as the final criteria
#' to determine whether or not the taxa should
#' be TRUE or FALSE. Just like with genefilter, a 
#' logical having length equal to nrow()/nspecies is returned, indicating which
#' should be kept. This output can be provided
#' directly to OTU trimming function, \code{\link{prune_species}}.
#'
#' @param X The object that needs trimming. Can be matrix, otuTable, or higher-
#' order phyloseq classes that contain an otuTable.
#'
#' @param flist An enclosure object, typically created with \code{\link{filterfunSample}}
#'
#' @param A An integer. The number of samples in which a taxa / species passed the filter
#' for it to be labeled TRUE in the output logical vector.
#'
#' @return A logical vector with names equal to species.names (or rownames, if matrix).
#'
#' @seealso \code{\link{prune_species}}
#' @keywords agglomerate OTU cluster tree
#'
#' @rdname genefilterSample-methods
#' @docType methods
#' @export
#'
#' @examples #
#' ## testOTU <- otuTable(matrix(sample(1:50, 25, replace=TRUE), 5, 5), speciesAreRows=FALSE)
#' ## f1  <- filterfunSample(topk(2))
#' ## wh1 <- genefilterSample(testOTU, f1, A=2)
#' ## wh2 <- c(T, T, T, F, F)
#' ## prune_species(wh1, testOTU)
#' ## prune_species(wh2, testOTU)
#' ## 
#' ## taxtab1 <- taxTab(matrix("abc", 5, 5))
#' ## prune_species(wh1, taxtab1)
#' ## prune_species(wh2, taxtab1)
setGeneric("genefilterSample", function(X, flist, A=1) standardGeneric("genefilterSample"))
#' @rdname genefilterSample-methods
#' @aliases genefilterSample,matrix-method
setMethod("genefilterSample", signature("matrix"), function(X, flist, A=1){
	TFmat = apply(X, 2, flist)
	apply(TFmat, 1, function(x, A){sum(x)>A}, A)
})
#' @rdname genefilterSample-methods
#' @aliases genefilterSample,otuTable-method
setMethod("genefilterSample", signature("otuTable"), function(X, flist, A=1){
	if( speciesAreRows(X) ){
		genefilterSample(   as(X, "matrix"), flist, A)
	} else {
		genefilterSample( t(as(X, "matrix")), flist, A)
	}
})
#' @rdname genefilterSample-methods
#' @aliases genefilterSample,phyloseqFather-method
setMethod("genefilterSample", signature("phyloseqFather"), function(X, flist, A=1){
	genefilterSample(otuTable(X), flist, A)
})
################################################################################
#' A sample-wise generic filter function, analogous to filterfun from the 
#' genefilter package. 
#'
#' @param ... A comma-separated list of functions.
#' 
#' @return An enclosure (function) that itself will return a logical vector, 
#'  according to the
#'  functions provided in the argument list, evaluated in order. The output of
#'  filterfunSample is appropriate for the `flist' argument to the 
#'  genefilterSample method.
#' 
#' @export
#' @seealso genefilter genefilterSample
#' @examples #
filterfunSample = function(...){
    flist <- list(...)
    if( length(flist) == 1 && is.list(flist[[1]])) { flist <- flist[[1]] }
    f = function(x){
    	# initialize fval (a logical vector)
    	fun  = flist[[1]]
    	fval = fun(x)
    	# check the remaining functions. Compare & logic, element-wise, each loop.
        for(fun in flist[-1]){
            fval = fval & fun(x)
		}
		return(fval)
	}
	class(f) <- "filterfun"
	return(f)
}
############################################################
#' The most abundant \code{k} taxa
#'
#' @param k An integer, indicating how many of the most abundant taxa
#'  should be kept.
#' @param na.rm A logical. Should \code{NA}s be removed. Default is \code{TRUE}.
#'
#' @return Returns a function (enclosure) that will return TRUE
#'  for each element in the most abundant k values.
#'
#' @export
topk = function(k, na.rm=TRUE){
    function(x){
		if(na.rm){x = x[!is.na(x)]}
		x >= sort(x, decreasing=TRUE)[k]
    }
}
############################################################
#' The most abundant \code{p} fraction of taxa
#'
#' @param p A numeric of length 1, indicating what fraction of the most abundant taxa
#'  should be kept.
#' @param na.rm A logical. Should \code{NA}s be removed. Default is \code{TRUE}.
#'
#' @return A function (enclosure), suitable for \code{\link{filterfunSample}},
#'  that will return \code{TRUE}
#'  for each element in the most abundant p fraction of taxa.
#'
#' @export
topp = function(p, na.rm=TRUE){
    function(x){
		if(na.rm){x = x[!is.na(x)]}
		x >= sort(x, decreasing=TRUE)[ceiling(length(x)*p)]
    }
}
################################################################################
#' The top f fraction of observations in a sample.
#'
#' As opposed to \code{\link{topp}}, which gives the
#' most abundant p fraction of observed taxa (richness, instead of cumulative
#' abundance. Said another way, topf ensures a certain
#' fraction of the total sequences are retained, while topp ensures
#' that a certain fraction of taxa/species/OTUs are retained.
#'
#' @param f Single numeric value between 0 and 1.
#' @param na.rm Logical. Should we remove NA values. Default \code{TRUE}.
#'
#' @return A function (enclosure), suitable for \code{\link{filterfunSample}},
#'  that will return \code{TRUE}
#'  for each element in the taxa comprising the most abundant f fraction of individuals.
#'
#' @export
#' @examples
#' t1 <- 1:10; names(t1)<-paste("t", 1:10, sep="")
#' topf(0.6)(t1)
topf <- function(f, na.rm=TRUE){
    function(x){
        if (na.rm){
            x = x[!is.na(x)]
        }
        y <- sort(x, decreasing = TRUE)
        y <- cumsum(y)/sum(x)
        return( (y <= f)[names(x)] )
    }
}
################################################################################
#' Set to FALSE any outlier species greater than f fractional abundance.
#'
#' This is for removing overly-abundant outlier taxa, not for trimming low-abundance
#' taxa.
#'
#' @param f Single numeric value between 0 and 1. The maximum fractional abundance
#'  value that a taxa will be allowed to have in a sample without being marked
#'  for trimming.
#'
#' @param na.rm Logical. Should we remove NA values. Default \code{TRUE}.
#'
#' @return A function (enclosure), suitable for \code{\link{filterfunSample}}.
#'
#' @export
#' @examples
#' t1 <- 1:10; names(t1)<-paste("t", 1:10, sep="")
#' rm_outlierf(0.15)(t1)
rm_outlierf <- function(f, na.rm=TRUE){
	function(x){
		if(na.rm){
			x = x[!is.na(x)]
		}
		y <- x / sum(x)
        return( y < f )
    }
}
################################################################################
# ################################################################################
# # A generic extension of the genefilter method from the genefilter package. 
# #
# # @param expr A phyloseq object that is or contains an otuTable.
# # @param flist A filterfun class object, a function, that is itself often
# # constructred from a list of functions. It takes a single argument, and
# # returns a logical vector with length equal to the number of rows in
# # \code{expr}.
# # 
# # @return A logical. See genefilter from the genefilter package.
# # 
# # @seealso genefilter genefilterSample
# # @keywords genefilter trim filter
# # @importFrom genefilter genefilter
# # @examples #
# setGeneric("genefilter")
# setMethod("genefilter", signature("phyloseqFather"),function(expr, flist){
	# genefilter(otuTable(expr), flist)
# })