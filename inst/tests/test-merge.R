# testthat tests don't do anything when successful.
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
test_that("merge_species() replaces disagreements in taxonomy with NA", {
	# Try a more difficult merge from a different subset
	GP20 <- prune_species(species.names(GlobalPatterns)[1:20], GlobalPatterns)
	
	# Arbitrary merge into taxa "951", NA in ranks after Phylum
	merge_these <- c("951", "586076", "141782", "30678", "30405")
	n5 <- merge_species(GP20, merge_these)
	# Test that none of the non-archetype taxa are left after the merge
	expect_that(all( !c("586076", "141782", "30678", "30405") %in% species.names(n5)), equals(TRUE))
	# Test that the archetype taxa remains
	expect_that( "951" %in% species.names(n5), equals(TRUE))
	# Test that the taxonomy is NA after Phylum
	n5_merged_taxonomy <- as(taxTab(n5), "matrix")["951", ]
	expect_that(any(is.na(n5_merged_taxonomy[1:2])), equals(FALSE))
	expect_that(all(is.na(n5_merged_taxonomy[3:7])), equals(TRUE))	
	
	# Test how well it works at a different level (say first or last ranks)
	merge_these <- c("1126", "31759")
	n6 <- merge_species(GP20, merge_these)
	# Test that the non-archetype taxa is gone
	expect_that( !"31759" %in% species.names(n6), equals(TRUE))
	# Test that the archetype taxa remains
	expect_that( "1126" %in% species.names(n6), equals(TRUE))
	# Test that the taxonomy is NA after Order
	n6_merged_taxonomy <- as(taxTab(n6), "matrix")["1126", ]
	expect_that( any(is.na(n6_merged_taxonomy[1:4])), equals(FALSE))
	expect_that( all(is.na(n6_merged_taxonomy[5:7])), equals(TRUE))	

	# Test that it works for differences at the first rank
	GP20f <- GP20
	taxTab(GP20f)[1, 1] <- "Bacteria"
	n7 <- merge_species(GP20f, species.names(GP20f)[1:2])
	# Should be all NA taxonomy
	expect_that( all(is.na(as(taxTab(n7), "matrix")[1, ])), equals(TRUE))

	# Test that it works for differences at the last rank
	# First, make the first taxa the same as "951"
	taxTab(GP20f)[1, ] <- taxTab(GP20f)["951", ]
	# Now change the last rank of this entry to something else
	taxTab(GP20f)[1, length(rank.names(GP20f))] <- "species_phyloseq_test"
	n8 <- merge_species(GP20f, c("951", species.names(GP20f)[1]))
	t951 <- as(taxTab(n8), "matrix")["951", ]	
	expect_that( sum(is.na(t951)), equals(1L))
	expect_that( is.na(t951[length(rank.names(n8))]), is_equivalent_to(TRUE))
	expect_that( t951[-7], is_identical_to(as(taxTab(GP20f), "matrix")["951", ][-7]))

	# Test that it works if the taxonomies completely agree
	GP20f <- GP20	
	# Make the first taxa the same as "951"
	taxTab(GP20f)[1, ] <- taxTab(GP20f)["951", ]
		
	merge_these <- c("549322", "951")
	n9   <- merge_species(GP20f, merge_these)
	n9t1 <- as(taxTab(n9), "matrix")["549322", ]
	# None should be NA
	expect_that(any(is.na(n9t1)), equals(FALSE))
	expect_that(length(n9t1), equals(7L))
	# Merge worked, "951" is gone.
	expect_that("951" %in% species.names(n9), equals(FALSE))	

})






