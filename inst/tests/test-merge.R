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

################################################################################
# taxglom
# Load data
data("GlobalPatterns")
GP.chl <- subset_species(GlobalPatterns, Phylum == "Chlamydiae")
test_that("the taxTab slot is identical whether taxglom()ed by itself or as component", {
	expect_that(taxglom(taxTab(GP.chl), "Family"), is_a("taxonomyTable"))
	expect_that(n1<-taxglom(GP.chl, "Family"), is_a("phyloseq"))
	expect_that(nspecies(n1), equals(4L))
	expect_that(
		taxglom(taxTab(GP.chl), "Family"),
		is_identical_to(taxTab(taxglom(GP.chl, "Family")))
	)
	expect_that(
		taxglom(taxTab(GP.chl), "Family", FALSE),
		is_identical_to(taxTab(n2<-taxglom(GP.chl, "Family", FALSE)))
	)
	expect_that(nspecies(n2), equals(5L))	
})
test_that("taxglom() handles clearly agglomeration to one taxa", {
	expect_that(n1 <- taxglom(GP.chl, "Phylum"), gives_warning())
	expect_that(n1, is_a("phyloseq"))
	expect_that(nspecies(n1), equals(1L))
	expect_that(access(n1, "tre"), is_a("NULL"))
})
################################################################################
# prune_species
# Use the GP.chl dataset from previous testing block
test_that("prune_species() handles clearly pruning to one taxa", {
	# throws warning, and NULL-tre
	expect_that(n1 <- prune_species(species.names(GP.chl)[1:1], GP.chl), gives_warning())
	expect_that(nspecies(n1), equals(1L))
	expect_that(n1, is_a("phyloseq"))
	expect_that(access(n1, "tre"), is_a("NULL"))
	expect_that(access(n1, "otuTable"), is_a("otuTable"))
})
test_that("prune_species() properly handles standard-cases", {
	# throws warning, and NULL-tre
	expect_that(n1 <- prune_species(species.names(GP.chl)[1:5], GP.chl), is_a("phyloseq"))
	expect_that(nspecies(n1), equals(5L))
	expect_that(access(n1, "tre"), is_a("phylo"))
	expect_that(access(n1, "otuTable"), is_a("otuTable"))
	expect_that(access(n1, "samData"), is_a("sampleData"))
	expect_that(access(n1, "taxTab"), is_a("taxonomyTable"))
	# Use logical vector, and get same answer
	L2 <- vector(length=nspecies(GP.chl))
	L2[1:5] <- TRUE
	expect_that(n2 <- prune_species(L2, GP.chl), is_a("phyloseq"))
	expect_that(n2, is_identical_to(n1))	
})
################################################################################
# merge_species
# Use the GP.chl dataset from previous testing block
test_that("merge_species() properly handles standard-cases", {
	expect_that(n1 <- merge_species(GP.chl, c("24341", "579085")), is_a("phyloseq"))
	expect_that(nspecies(n1), equals(20L))
	# The first name is kept, others removed
	expect_that("579085" %in% species.names(n1), equals(FALSE))
	expect_that("24341" %in% species.names(n1),  equals(TRUE))
	# Try a 3-element merge
	expect_that(n2 <- merge_species(GP.chl, c("579085", "24341", "547579")), is_a("phyloseq"))
	expect_that(nspecies(n2), equals(19L))
	# The first name is kept, others removed
	expect_that("579085" %in% species.names(n2), equals(TRUE))
	expect_that("24341"  %in% species.names(n2), equals(FALSE))
	expect_that("547579" %in% species.names(n2), equals(FALSE))	
	# Try again, but specify the retained OTU name as the 3rd one
	expect_that(n3 <- merge_species(GP.chl, c("579085", "24341", "547579"), "547579"), is_a("phyloseq"))
	# "547579" is kept, others removed
	expect_that("579085" %in% species.names(n3), equals(FALSE))
	expect_that("24341"  %in% species.names(n3), equals(FALSE))
	expect_that("547579" %in% species.names(n3), equals(TRUE))
	# Check that the remaining OTU has the sum of the values merged
	expect_that(getSamples(n3, "547579"), is_identical_to(
			getSamples(GP.chl, "24341") + getSamples(GP.chl, "547579") + getSamples(GP.chl, "579085")
		)
	)
})
# Need to check the taxonomyTable results from these merges...
# taxTab(n3)["547579", ]
# taxTab(GP.chl)["24341", ]





