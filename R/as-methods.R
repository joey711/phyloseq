################################################################################
# coercion methods
################################################################################
setAs("phyloseq", "matrix", function(from){
	from@.Data
})
setAs("phyloseq", "otu_table", function(from){
	otu_table(from)	
})
setAs("phyloseq", "otu_table", function(from){
	otu_table(from)	
})
################################################################################
setAs("data.frame", "sample_data", function(from){
	new("sample_data", from)
})
setAs("sample_data", "data.frame", function(from){
	data.frame(from)
})
setAs("phyloseq", "sample_data", function(from){
	sample_data(from)	
})
################################################################################
setAs("taxonomyTable", "matrix", function(from){
	from@.Data
})
setAs("phyloseq", "taxonomyTable", function(from){
	tax_table(from)	
})
################################################################################
setAs("phyloseq", "phylo", function(from){
	phy_tree(from)
})
################################################################################
