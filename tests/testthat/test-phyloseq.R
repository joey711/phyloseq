################################################################################
# Use testthat to test phyloseq constructor and other class internals.
################################################################################
library("phyloseq"); library("testthat")
# # # # TESTS!
set.seed(8888)

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
	expect_identical(esophagus1, esophagus)
	# Now try to "rebuild" using the quote-containing 
	esophagus2 = phyloseq(tree, otu_table(esophagus)) 
	expect_identical(esophagus2, esophagus)
	
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
	expect_identical(GP1, GP)
	# Now try to "rebuild" using the quote-containing 
	GP2 = phyloseq(tree, otu_table(GP), tax_table(GP), sample_data(GP)) 
	expect_identical(GP2, GP)
})
################################################################################
# More constructor-related tests needed here...

################################################################################
# - Test that re-assigning taxa_names and sample_names works.
# - Use this to test that intersect_samples and intersect_taxa works.
################################################################################
data("GlobalPatterns")
e1 = prune_taxa(taxa_names(GlobalPatterns)[1:25], GlobalPatterns)

test_that("taxa_names(x)<- and sample_names(x)<- behaves as expected", {
	# taxa_names<-
	new_taxa_names = paste("OTU-", taxa_names(e1), sep="")
	taxa_names(e1) = new_taxa_names
	expect_true(identical(taxa_names(e1), new_taxa_names))
	expect_true(identical(taxa_names(phy_tree(e1)), new_taxa_names))
	expect_true(identical(taxa_names(otu_table(e1)), new_taxa_names))
	expect_true(identical(taxa_names(tax_table(e1)), new_taxa_names))
# 	expect_true(identical(taxa_names(refseq(e1)), new_taxa_names))
	
	# sample_names<-
	new_sample_names = paste("Sa-", sample_names(e1), sep="")
	sample_names(e1) = new_sample_names
	expect_true(identical(sample_names(e1), new_sample_names))
	expect_true(identical(sample_names(sample_data(e1)), new_sample_names))
	expect_true(identical(sample_names(otu_table(e1)), new_sample_names))
})
test_that("Test intersect_*() and prune_*() methods behave as expected", {
	e0 = e1
	# taxa_names<-
	## We assign new names to just one component, being sneaky and using the
	## not-recommended direct replacement with @slotname
	## This should work, but users should not do it in normal circumstances
	new_taxa_names = taxa_names(e1)
	nunchained = 5L
	i = sample(ntaxa(e1), ntaxa(e1)-nunchained, replace=FALSE)
	new_taxa_names[i] = paste("OTU-", taxa_names(e1)[i], sep="")
	taxa_names(e1@tax_table) = new_taxa_names
	expect_false(identical(taxa_names(e1), new_taxa_names))
	expect_true(identical(taxa_names(tax_table(e1)), new_taxa_names))	
	expect_false(identical(taxa_names(otu_table(e1)), new_taxa_names))	
	expect_false(identical(taxa_names(phy_tree(e1)), new_taxa_names))	
# 	expect_false(identical(taxa_names(refseq(e1)), new_taxa_names))	
	## Okay so that worked. Now we test if the intersection functions behave
	expect_true(identical(length(phyloseq:::intersect_taxa(e1)), nunchained))
	e2 = prune_taxa(phyloseq:::intersect_taxa(e1), e1)
	expect_true(identical(ntaxa(e2), nunchained))
	expect_true(setequal(taxa_names(e2), taxa_names(e1)[-i]))
	
	# sample_names<-
	e1 = e0
	## We assign new names to just one component, being sneaky and using the
	## not-recommended direct replacement with @slotname
	## This should work, but users should not do it in normal circumstances
	new_sample_names = sample_names(e1)
	nunchained = 5L
	i = sample(nsamples(e1), nsamples(e1)-nunchained, replace=FALSE)
	new_sample_names[i] = paste("Sa-", sample_names(e1)[i], sep="")
	sample_names(e1@sam_data) = new_sample_names
	expect_false(identical(sample_names(e1), new_sample_names))
	expect_true(identical(sample_names(sample_data(e1)), new_sample_names))	
	expect_false(identical(sample_names(otu_table(e1)), new_sample_names))	
	
	## Okay so that worked. Now we test if the intersection functions behave
	expect_true(identical(length(phyloseq:::intersect_samples(e1)), nunchained))
	e2 = prune_samples(phyloseq:::intersect_samples(e1), e1)
	expect_true(identical(nsamples(e2), nunchained))
	expect_true(setequal(sample_names(e2), sample_names(e1)[-i]))
	
})

test_that("Test ordering", {
	OTU = otu_table(e1)
	tree = phy_tree(e1)
	expect_true(identical(taxa_names(OTU), taxa_names(tree)))
	reotaxnames = sample(taxa_names(tree), ntaxa(tree), FALSE)
	expect_false(identical(taxa_names(OTU), reotaxnames))
	# scramble order of taxa_names in tree by random arbitrary assignment
	taxa_names(tree) <- reotaxnames
	expect_false(identical(taxa_names(OTU), taxa_names(tree)))
	# implicitly re-order in constructor
	e3 = phyloseq(OTU, tree)
	expect_true(identical(taxa_names(e3), taxa_names(tree)))
	expect_true(identical(taxa_names(otu_table(e3)), taxa_names(phy_tree(e3))))
	expect_true(identical(taxa_names(otu_table(e3)), taxa_names(phy_tree(e3))))
	
	# Glad that worked. Now let's mess up sample indices in one component, and OTU indices in another
	# then fix explicitly with index_reorder(e4, "both")
	e4 = e1
	reosamplenames = sample(sample_names(e1), nsamples(e1), FALSE)
	sample_names(e4@sam_data) <- reosamplenames
	taxa_names(e4@tax_table) <- reotaxnames
	
	expect_false(identical(taxa_names(otu_table(e4)), taxa_names(tax_table(e4))))
	expect_false(identical(sample_names(otu_table(e4)), sample_names(sample_data(e4))))
	e4 = phyloseq:::index_reorder(e4, "both")
	expect_true(identical(taxa_names(otu_table(e4)), taxa_names(tax_table(e4))))
	expect_true(identical(sample_names(otu_table(e4)), sample_names(sample_data(e4))))	
})

################################################################################