################################################################################
# Use testthat to test phyloseq transformation functions/methods
################################################################################
library("phyloseq"); library("testthat")
# # # # TESTS!
################################################################################
# rarefy_even_depth
################################################################################
data("GlobalPatterns")
set.seed(711) # The random seed for randomly selecting subset of OTUs
randoOTUs = sample(taxa_names(GlobalPatterns), 100, FALSE)
GP100 = prune_taxa(randoOTUs, GlobalPatterns)
min_lib = 1000
# The default rng seed is being implied in this call (also 711)
rGP  = suppressMessages(rarefy_even_depth(GP100, sample.size=min_lib, rngseed=FALSE))
rGPr = suppressMessages(rarefy_even_depth(GP100, sample.size=min_lib, rngseed=FALSE, replace=FALSE))
################################################################################
# Test that specific OTUs and samples were removed
################################################################################
test_that("Test that empty OTUs and samples were automatically pruned", {
	rmOTU = setdiff(taxa_names(GP100), taxa_names(rGP))
	expect_equal(length(rmOTU), 20L)
	expect_equal(rmOTU[1:5], c("534601", "408325", "325564", "8112", "571917"))
	expect_true(taxa_names(GP100)[taxa_sums(GP100) <= 0] %in% rmOTU)
	expect_true(all(taxa_sums(rGP) > 0))
	rmsam = setdiff(sample_names(GP100), sample_names(rGP))
	expect_equal(length(rmsam), 12L)
	expect_equal(rmsam[1:5], c("M11Fcsw", "M31Tong", "M11Tong", "NP2", "TRRsed1"))
	expect_true(all(sample_sums(rGP) > 0))
})
################################################################################
# Test specific values. Should be reproducible, and you set the seed.
################################################################################
test_that("Test values", {
  # with replacement values
	expect_equal(as(otu_table(rGP)[1, 3:10], "vector"), rep(0, 8))
	expect_equal(as(otu_table(rGP)[2, 1:10], "vector"), c(rep(0, 9), 2))
	expect_equal(as(otu_table(rGP)[3, 8:12], "vector"), c(892, 956, 56, 10, 25))
	expect_equal(as(otu_table(rGP)[70:78, 4], "vector"),
							 c(710, 2, 0, 2, 0, 8, 154, 2, 0))
	# without replacement values
	expect_equal(as(otu_table(rGPr)[1, 3:10], "vector"), c(rep(0, 7), 1))
	expect_equal(as(otu_table(rGPr)[2, 1:10], "vector"), 
               c(rep(0, 5), 4, 0, 877, 960, 55))
	expect_equal(as(otu_table(rGPr)[3, 8:12], "vector"), 
               c(10, 34, 2, 0, 2))
	expect_equal(as(otu_table(rGPr)[70:78, 4], "vector"),
	             c(0, 706, 1, 0, 2, 0, 5, 173, 1))  
})
################################################################################
# Include tests from the rarefy-without-replacement results, used by many.
#################################################################################
test_that("Test library sizes are all the same set value", {
  expect_true(all(sample_sums(rGP )==min_lib))
  expect_true(all(sample_sums(rGPr)==min_lib))
})
test_that("The same samples should have been cut in each results", {
  expect_equal(nsamples(rGP), 14)
  expect_true(setequal(sample_names(rGP), sample_names(rGPr)))
})
################################################################################