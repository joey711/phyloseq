################################################################################
# coercion methods
################################################################################
setAs("phyloseq", "matrix", function(from){
	from@.Data
})
setAs("phyloseq", "otuTable", function(from){
	otuTable(from)	
})
setAs("phyloseq", "otuTable", function(from){
	otuTable(from)	
})
################################################################################
setAs("data.frame", "sampleData", function(from){
	new("sampleData", from)
})
setAs("sampleData", "data.frame", function(from){
	data.frame(from)
})
setAs("phyloseq", "sampleData", function(from){
	sampleData(from)	
})
################################################################################
setAs("taxonomyTable", "matrix", function(from){
	from@.Data
})
setAs("phyloseq", "taxonomyTable", function(from){
	taxTab(from)	
})
################################################################################
setAs("phyloseq", "phylo", function(from){
	tre(from)
})
################################################################################