################################################################################
# coercion methods
# Google how to document setAs statements using roxygen.
################################################################################
#' Coerce an otuTable object to a matrix.
#'
#' @name as
#' @aliases as,otuTable,matrix-method
#' @docType methods
#' @rdname as-methods
setAs("otuTable", "matrix", function(from){
	from@.Data
})
setAs("otuTree", "otuTable", function(from){
	otuTable(from)	
})
setAs("otuSam", "otuTable", function(from){
	otuTable(from)	
})
################################################################################
setAs("sampleMap", "data.frame", function(from){
	new("sampleMap", from)
})
setAs("data.frame", "sampleMap", function(from){
	data.frame(from)
})
setAs("otuSam", "sampleMap", function(from){
	sampleMap(from)	
})
################################################################################
setAs("taxonomyTable", "matrix", function(from){
	from@.Data
})
setAs("otuSamTax", "taxonomyTable", function(from){
	taxTab(from)	
})
################################################################################
