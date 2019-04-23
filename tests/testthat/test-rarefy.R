################################################################################
# Use testthat to test phyloseq transformation functions/methods
################################################################################
library("phyloseq"); library("magrittr"); library("testthat")
# # # # TESTS!
################################################################################
# rarefy_even_depth
################################################################################
data("GlobalPatterns")
# The random seed for randomly selecting subset of OTUs
rngseed = 711 
keepTaxa <- 
  GlobalPatterns %>% taxa_sums() %>% sort %>% tail(100) %>% names() %>% 
  # Append some OTUs that will be cut...
  # c(., (GlobalPatterns %>% taxa_sums() %>% sort %>% names() %>% .[101:length(.)]) %>% sample(size = 100))
  c(., c("12589", "10444", "9286", "374370", "63062", "158132", "324145", "180450", "178513", "542714"))
GP100 = prune_taxa(keepTaxa, GlobalPatterns)
min_lib = 4000
# The default rng seed is being implied in this call
rGP  = rarefy_even_depth(GP100, sample.size=min_lib, rngseed=rngseed)
rGPr = rarefy_even_depth(GP100, sample.size=min_lib, rngseed=rngseed + 1L, replace=FALSE)
################################################################################
# Test that specific OTUs and samples were removed
################################################################################
test_that("Test that empty OTUs and samples were automatically pruned", {
	rmOTU = setdiff(taxa_names(GP100), taxa_names(rGP))
	setdiff(taxa_names(rGP), taxa_names(rGPr))
	expect_equal(length(rmOTU), 3L)
	expect_equal(rmOTU, c("10444", "9286", "542714"))
	expect_true(all(taxa_sums(rGP) > 0))
	rmsam = setdiff(sample_names(GP100), sample_names(rGP))
	expect_equal(length(rmsam), 2L)
	expect_equal(rmsam, c("CC1", "SV1"))
	expect_true(all(sample_sums(rGP) > 0))
})
################################################################################
# Include tests from the rarefy-without-replacement results, used by many.
#################################################################################
test_that("Test library sizes are all the same set value", {
  expect_true(all(sample_sums(rGP )==min_lib))
  expect_true(all(sample_sums(rGPr)==min_lib))
})
test_that("The same samples should have been cut in each results", {
  expect_equal(nsamples(rGP), 24)
  expect_true(setequal(sample_names(rGP), sample_names(rGPr)))
})
################################################################################
