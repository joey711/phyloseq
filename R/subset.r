################################################################################
# subsetting functions
# Want subset of our special tables to still be those classes after subset
# Without these, the default converts to the base object (e.g. matrix or data.frame)
################################################################################
setMethod("[", "otuTable", function(x,i,j,...){
	newx <- callNextMethod(x@.Data,i,j,drop=FALSE,...)
	new("otuTable", x@speciesAreRows, newx)
})
################################################################################
setMethod("[", "sampleMap", function(x,i,j,...){
	new("sampleMap", callNextMethod(data.frame(x),i,j,drop=FALSE,...))
})
################################################################################
setMethod("[", "taxonomyTable", function(x,i,j,...){
	new("taxonomyTable", callNextMethod(x@.Data,i,j,drop=FALSE,...))
})
################################################################################
################################################################################
