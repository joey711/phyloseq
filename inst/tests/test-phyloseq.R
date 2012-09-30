################################################################################
# Use testthat to test phyloseq constructor and other class internals.
################################################################################
library("phyloseq"); library("testthat")
# # # # TESTS!

################################################################################
test_that("phyloseq: Building a phyloseq-object when tree contains extra quotes, still works.", {
	data("esophagus")
	tree = phy_tree(esophagus)
	# Add extra quotes surrounding each OTU name in the tree
	tree$tip.label = paste("\"", taxa_names(tree), "\"", sep="")
	# Try to add the tree back to esophagus, replacing the original
	# (Should work with message.)
	esophagus1 = esophagus
	phy_tree(esophagus) = tree
	expect_that(esophagus1, is_identical_to(esophagus))
	# Now try to "rebuild" using the quote-containing 
	esophagus2 = phyloseq(tree, otu_table(esophagus)) 
	expect_that(esophagus2, is_identical_to(esophagus))
	# Check for message:
	expect_that(phyloseq(tree, otu_table(esophagus)), shows_message())
	
	# Try with a dataset with complete-set of components, Global Patterns
	data("GlobalPatterns")
	# Use a subset, because checking identicality of two large, complicated objects takes time.
	minsum = sort(taxa_sums(GlobalPatterns), TRUE)[20]
	GP = prune_taxa(taxa_sums(GlobalPatterns) >= minsum, GlobalPatterns)
	tree = phy_tree(GP)
	# Add extra quotes surrounding each OTU name in the tree
	tree$tip.label = paste("\"", taxa_names(tree), "\"", sep="")
	# Try to add the tree back to GP, replacing the original
	# (Should work with message.)
	GP1 = GP
	phy_tree(GP) = tree
	expect_that(GP1, is_identical_to(GP))
	# Now try to "rebuild" using the quote-containing 
	GP2 = phyloseq(tree, otu_table(GP), tax_table(GP), sample_data(GP)) 
	expect_that(GP2, is_identical_to(GP))
	# Check for message:
	expect_that(phyloseq(tree, otu_table(GP)), shows_message())	
})
################################################################################
# More constructor-related tests needed here...

