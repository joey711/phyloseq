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
setAs("data.frame", "sampleMap", function(from){
	new("sampleMap", from)
})
setAs("sampleMap", "data.frame", function(from){
	data.frame(from)
})
setAs("phyloseq", "sampleMap", function(from){
	sampleMap(from)	
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