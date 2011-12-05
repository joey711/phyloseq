################################################################################
#' Merge samples based on a sample variable or factor.
#'
#' The purpose of this method is to merge/agglomerate the sample indices of a 
#' phyloseq object according to a categorical variable contained in a sampleMap
#' or a provided factor.
#' 
#' NOTE: (\code{\link[phylobase]{phylo4}}) trees and \code{\link{taxonomyTable-class}}
#' are not modified by this function, but returned in the output object as-is. 
#'
#' @usage merge_samples(x, group, fun=mean) 
#'
#' @param x (Required). An instance of a phyloseq class that has sample indices. This includes 
#'  \code{\link{sampleMap-class}}, \code{\link{otuTable-class}}, and \code{\link{phyloseq-class}}. 
#'
#' @param group (Required). Either the a single character string matching a variable name in
#'  the corresponding sampleMap of \code{x}, or a factor with the same length as
#'  the number of samples in \code{x}.
#'
#' @param fun (Optional). The function that will be used to merge the values that
#'  correspond to the same group for each variable. It must take a numeric vector
#'  as first argument and return a single value. Default is \code{\link[base]{mean}}.
#'  Note that this is (currently) ignored for the otuTable, where the equivalent
#'  function is \code{\link[base]{sum}}, but evaluated via \code{\link[base]{rowsum}}
#'  for efficiency.
#'
#' @return A phyloseq object that has had its sample indices merged according to
#'  the factor indicated by the \code{group} argument. The output class
#'  matches \code{x}.  
#'
#' @seealso \code{\link{mergespecies}}, code{\link{merge_phyloseq}}
#'
#' @rdname merge_samples-methods
#' @docType methods
#' @export
#'
#' @examples #
#' # # # data(ex1)
#' # # # t1 <- merge_samples(sampleMap(ex1), "Gender")
#' # # # t4 <- merge_samples(ex1, "Gender")
#' # # # identical(t1, sampleMap(t4))
#' # # # t0 <- merge_samples(ex1, "Diet")
#' # # # grp <- as(data.frame(sampleMap(ex1))[, "Diet"], "vector")
#' # # # t4 <- merge_samples(ex1, grp)
#' # # # identical(t0, t4)
#' # # # t1 <- merge_samples(otuTable(ex1), grp)
#' # # # t2 <- merge_samples(otuTable(ex1), factor(grp))
#' # # # t3 <- merge_samples(ex1, "Diet")
#' # # # identical(t1, t2)
#' # # # identical(t1, t3)
#' # # # identical(t1, otuTable(t3))
setGeneric("merge_samples", function(x, group, fun=mean) standardGeneric("merge_samples"))
################################################################################
#' @aliases merge_samples,sampleMap-method
#' @rdname merge_samples-methods
setMethod("merge_samples", signature("sampleMap"), function(x, group, fun=mean){
	x1    <- data.frame(x)

	# Check class of group and modify if "character"
	if( class(group)=="character" & length(group)==1 ){
		if( !group %in% colnames(x) ){stop("group not found among sample variable names.")}
		group <- x1[, group]
	}
	if( class(group)!="factor" ){
		# attempt to coerce to factor
		group <- factor(group)
	}

	# Remove any non-coercable columns.
	# Philosophy is to keep as much as possible. If it is coercable at all, keep.
	# Coerce all columns to numeric matrix
	coercable    <- sapply(x1, canCoerce, "numeric")
	x2           <- sapply(x1[, coercable], as, "numeric")
	rownames(x2) <- rownames(x1)	
	
	# Perform the aggregation.
	outdf <- aggregate(x2, list(group), fun)
	# get rownames from the "group" column (always first)
	# rownames(outdf) <- as.character(outdf[, 1])
	rownames(outdf) <- levels(group)
	# "pop" the first column
	outdf <- outdf[, -1, drop=FALSE]

	return( sampleMap(outdf) )
})
################################################################################
#' @aliases merge_samples,otuTable-method
#' @rdname merge_samples-methods
setMethod("merge_samples", signature("otuTable"), function(x, group){
	# needs to be in sample-by-species orientation
	if( speciesAreRows(x) ){ x <- t(x) }
	# coerce to matrix, x2
	x2 <- as(x, "matrix")
	
	# # # #aggregate(x2, list(group), fun)
	out <- rowsum(x2, group)
	
	# convert back to otuTable, and return
	return( otuTable(out, speciesAreRows=FALSE) )
})
################################################################################
#' @aliases merge_samples,phyloseq-method
#' @rdname merge_samples-methods
setMethod("merge_samples", signature("phyloseq"), function(x, group, fun=mean){

	# Check if phyloseq object has a sampleMap
	if( !is.null(access(x, "sampleMap")) ){
		# Check class of group and modify if single "character" (column name)
		if( class(group)=="character" & length(group)==1 ){
			x1 <- data.frame(sampleMap(x))		
			if( !group %in% colnames(x1) ){stop("group not found among sample variable names.")}
			group <- x1[, group]
		}
		if( class(group)!="factor" ){
			# attempt to coerce to factor
			group <- factor(group)
		}
		newSM <- merge_samples(sampleMap(x), group, fun)
		newOT <- merge_samples(otuTable(x), group)
		phyloseqList <- list(newOT, newSM)
		
	# Else, the only relevant object to "merge_samples" is the otuTable
	} else {
		if( class(group)!="factor" ){ group <- factor(group) }
		phyloseqList <- list( newOT=merge_samples(otuTable(x), group) )
	}
	
	### Add to build-call-list the remaining components, if present in x.
	### NULL is returned by accessor if object lacks requested component/slot.
	### Order of objects in list doesn't matter for phyloseq.
	### The list should not be named.
	if( !is.null(access(x, "taxTab")) ){ phyloseqList <- c(phyloseqList, list(taxTab(x))) }
	if( !is.null(access(x, "tre"))    ){ phyloseqList <- c(phyloseqList, list(tre(x))) }
	
	return( do.call("phyloseq", phyloseqList) )
})
################################################################################