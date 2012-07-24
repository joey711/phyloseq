# load libraries
library("phyloseq"); library("testthat")
# # # # TESTS!

# Load GP dataset
data("GlobalPatterns")
GP <- GlobalPatterns
keepNames <- sample.names(GP)[5:7]

test_that("Classes of pruned phyloseq and its components are as expected", {
	GP3     <- prune_samples(keepNames, GP)
	expect_that(nsamples(GP3), is_identical_to(3L))
	expect_that(GP3, is_a("phyloseq"))
	expect_that(access(GP3, "samData"),  is_a("sampleData"))
	expect_that(access(GP3, "otuTable"), is_a("otuTable"))	
	expect_that(access(GP3, "tre"),      is_a("phylo"))	
	expect_that(access(GP3, "taxTab"),   is_a("taxonomyTable"))		
	# Now try on instance without sample data (empty slot)
	GPnoSD  <- phyloseq(otuTable(GP), taxTab(GP))
	GP3noSD <- prune_samples(keepNames, GPnoSD)
	expect_that(nsamples(GP3noSD), is_identical_to(3L))	
	expect_that(access(GP3noSD, "otuTable"), is_a("otuTable"))	
	expect_that(access(GP3noSD, "samData"),  is_a("NULL"))
	expect_that(access(GP3noSD, "tre"),      is_a("NULL"))	
	expect_that(access(GP3noSD, "taxTab"),   is_a("taxonomyTable"))	
})

test_that("prune_samples works on sampleData-only and otuTable-only data", {
	GPotu <- prune_samples(keepNames, access(GP, "otuTable", TRUE))
	GPsd  <- prune_samples(keepNames, access(GP, "samData", TRUE))	
	expect_that(nsamples(GPotu), is_identical_to(3L))	
	expect_that(nsamples(GPsd), is_identical_to(3L))			
	expect_that(GPotu, is_a("otuTable"))
	expect_that(GPsd, is_a("sampleData"))
	expect_that(dim(GPotu), is_identical_to(c(19216L, 3L)))
	expect_that(dim(GPsd), is_identical_to(c(3L, 7L)))			
})

################################################################################