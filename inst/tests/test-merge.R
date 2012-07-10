# The following code should work, and return TRUE at the end
library("phyloseq")
library("testthat")


# # # Tests!
data(GlobalPatterns)
# GP <- prune_species(speciesSums(GlobalPatterns)>0, GlobalPatterns)
GP  <- GlobalPatterns
mGP <- merge_samples(GlobalPatterns, "SampleType")

test_that("Classes of merged phyloseq objects are as expected", {
	expect_that(merge_samples(otuTable(GP), getVariable(GP, "SampleType")), is_a("otuTable"))
	expect_that(merge_samples(sampleData(GP), "SampleType"), is_a("sampleData"))
	expect_that(mGP, is_a("phyloseq"))
})

test_that("Same samData result for separate and combined merge in merge_samples", {
	expect_that(
		merge_samples(sampleData(GP), "SampleType"),
		is_identical_to(sampleData(mGP))
	)
})

test_that("Same otuTable result for separate and combined merge in merge_samples", {
	expect_that(
		merge_samples(otuTable(GP), getVariable(GP, "SampleType")),
		is_identical_to(otuTable(mGP))
	)
})
