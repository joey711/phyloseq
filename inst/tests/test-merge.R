# The following code should work, and return TRUE at the end
library("phyloseq")
library("testthat")


# # # Tests!

################################################################################
# merge_samples
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

################################################################################
# merge_phyloseq
test_that("merge_phyloseq: Break apart GP based on human-association, then merge back together.", {
	data(GlobalPatterns)
	GP  <- prune_species(species.names(GlobalPatterns)[1:100], GlobalPatterns)
	sampleData(GP)$human <- factor(getVariable(GP, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue"))
	h1 <- subset_samples(GP, human=="TRUE")
	h0 <- subset_samples(GP, human=="FALSE")
	GP1 <- merge_phyloseq(h0, h1)

	# The species order is fixed to the tree, so should be the same between the original and merged
	expect_that(species.names(GP), is_identical_to(species.names(GP1)))
	expect_that(tre(h1), is_identical_to(tre(h0)))

	# However, the sample order has been shuffled by the split/merge. 
	# Fix the sample order by re-ordering the otuTable, and reassigning
	sa.order <- sample.names(GP)
	sa.order <- sa.order[sa.order %in% sample.names(GP1)]
	otuTable(GP1) <- otuTable(GP1)[, sa.order]

	# Should be fixed now. Full object and components now identical
	expect_that(GP1, equals(GP)) 
	expect_that(GP1, is_identical_to(GP))
	expect_that(otuTable(GP1), is_identical_to(otuTable(GP)))
	expect_that(sampleData(GP1), is_identical_to(sampleData(GP)))
	expect_that(taxTab(GP1), is_identical_to(taxTab(GP)))
	expect_that(tre(GP1), is_identical_to(tre(GP)))
})

