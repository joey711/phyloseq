################################################################################
# coercion methods
################################################################################
setAs("otuTable", "matrix",function(from){
	from@.Data
})
setAs("otuTree", "otuTable", function(from){
	otuTable(from)	
})
setAs("phyloseq", "otuTable", function(from){
	otuTable(from)	
})
################################################################################
setAs("sampleMap", "data.frame", function(from){
	new("sampleMap", from)
})
setAs("data.frame", "sampleMap", function(from){
	data.frame(from)
})
setAs("phyloseq", "sampleMap", function(from){
	sampleMap(from)	
})
################################################################################
setAs("taxonomyTable", "matrix", function(from){
	from@.Data
})
setAs("phyloseqTax", "taxonomyTable", function(from){
	taxTab(from)	
})
################################################################################
