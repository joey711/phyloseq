# testthat tests don't do anything when successful.
library("phyloseq")
library("testthat")


# # # Tests!

################################################################################
# merge_samples
data(GlobalPatterns)
# GP <- prune_taxa(taxa_sums(GlobalPatterns)>0, GlobalPatterns)
GP  <- GlobalPatterns
mGP <- merge_samples(GlobalPatterns, "SampleType")

test_that("Classes of merged phyloseq objects are as expected", {
	expect_that(merge_samples(otu_table(GP), get_variable(GP, "SampleType")), is_a("otu_table"))
	expect_that(merge_samples(sample_data(GP), "SampleType"), is_a("sample_data"))
	expect_that(mGP, is_a("phyloseq"))
})

test_that("Same sam_data result for separate and combined merge in merge_samples", {
	expect_that(
		merge_samples(sample_data(GP), "SampleType"),
		is_identical_to(sample_data(mGP))
	)
})

test_that("Same otu_table result for separate and combined merge in merge_samples", {
	expect_that(
		merge_samples(otu_table(GP), get_variable(GP, "SampleType")),
		is_identical_to(otu_table(mGP))
	)
})

test_that("Sample Names of merged object now same set as merging factor levels", {
	sampleTypes = levels(data.frame(sample_data(GP))$SampleType)
	expect_that(setdiff(sampleTypes, sample_names(mGP)), is_identical_to(character()))
})

test_that("Counts from merged-samples are summed...", {
	OTUnames10 = names(sort(taxa_sums(GP), TRUE)[1:10])
	GP10  = prune_taxa(OTUnames10,  GP)
	mGP10 = prune_taxa(OTUnames10, mGP)
	# Loop to check the correct summation has occured for all OTUs.
	for( i in OTUnames10 ){
		isum = as(tapply(get_sample(GP10, i), get_variable(GP10, "SampleType"), sum), "numeric")
		expect_that(isum, is_equivalent_to(get_sample(mGP10, i)))
	}
})

################################################################################
# merge_phyloseq
test_that("merge_phyloseq: Break apart GP based on human-association, then merge back together.", {
	data(GlobalPatterns)
	GP  <- prune_taxa(taxa_names(GlobalPatterns)[1:100], GlobalPatterns)
	sample_data(GP)$human <- factor(get_variable(GP, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue"))
	h1 <- subset_samples(GP, human=="TRUE")
	h0 <- subset_samples(GP, human=="FALSE")
	GP1 <- merge_phyloseq(h0, h1)

	# The species order is fixed to the tree, so should be the same between the original and merged
	expect_that(taxa_names(GP), is_identical_to(taxa_names(GP1)))
	expect_that(phy_tree(h1), is_identical_to(phy_tree(h0)))

	# However, the sample order has been shuffled by the split/merge. 
	# Fix the sample order by re-ordering the otu_table, and reassigning
	sa.order <- sample_names(GP)
	sa.order <- sa.order[sa.order %in% sample_names(GP1)]
	otu_table(GP1) <- otu_table(GP1)[, sa.order]

	# Should be fixed now. Full object and components now identical
	expect_that(GP1, equals(GP)) 
	expect_that(GP1, is_identical_to(GP))
	expect_that(otu_table(GP1), is_identical_to(otu_table(GP)))
	expect_that(sample_data(GP1), is_identical_to(sample_data(GP)))
	expect_that(tax_table(GP1), is_identical_to(tax_table(GP)))
	expect_that(phy_tree(GP1), is_identical_to(phy_tree(GP)))
})

################################################################################
# tax_glom
# Load data
data("GlobalPatterns")
GP.chl <- subset_taxa(GlobalPatterns, Phylum == "Chlamydiae")
test_that("the tax_table slot is identical whether tax_glom()ed by itself or as component", {
	expect_that(tax_glom(tax_table(GP.chl), "Family"), is_a("taxonomyTable"))
	expect_that(n1<-tax_glom(GP.chl, "Family"), is_a("phyloseq"))
	expect_that(ntaxa(n1), equals(4L))
	expect_that(
		tax_glom(tax_table(GP.chl), "Family"),
		is_identical_to(tax_table(tax_glom(GP.chl, "Family")))
	)
	expect_that(
		tax_glom(tax_table(GP.chl), "Family", FALSE),
		is_identical_to(tax_table(n2<-tax_glom(GP.chl, "Family", FALSE)))
	)
	expect_that(ntaxa(n2), equals(5L))	
})
test_that("tax_glom() handles clearly agglomeration to one taxa", {
	expect_that(n1 <- tax_glom(GP.chl, "Phylum"), gives_warning())
	expect_that(n1, is_a("phyloseq"))
	expect_that(ntaxa(n1), equals(1L))
	expect_that(access(n1, "phy_tree"), is_a("NULL"))
})
################################################################################
# prune_taxa
# Use the GP.chl dataset from previous testing block
test_that("prune_taxa() handles clearly pruning to one taxa", {
	# throws warning, and NULL-tre
	expect_that(n1 <- prune_taxa(taxa_names(GP.chl)[1:1], GP.chl), gives_warning())
	expect_that(ntaxa(n1), equals(1L))
	expect_that(n1, is_a("phyloseq"))
	expect_that(access(n1, "phy_tree"), is_a("NULL"))
	expect_that(access(n1, "otu_table"), is_a("otu_table"))
})
test_that("prune_taxa() properly handles standard-cases", {
	# throws warning, and NULL-tre
	expect_that(n1 <- prune_taxa(taxa_names(GP.chl)[1:5], GP.chl), is_a("phyloseq"))
	expect_that(ntaxa(n1), equals(5L))
	expect_that(access(n1, "phy_tree"), is_a("phylo"))
	expect_that(access(n1, "otu_table"), is_a("otu_table"))
	expect_that(access(n1, "sam_data"), is_a("sample_data"))
	expect_that(access(n1, "tax_table"), is_a("taxonomyTable"))
	# Use logical vector, and get same answer
	L2 <- vector(length=ntaxa(GP.chl))
	L2[1:5] <- TRUE
	expect_that(n2 <- prune_taxa(L2, GP.chl), is_a("phyloseq"))
	expect_that(n2, is_identical_to(n1))	
})
################################################################################
# merge_taxa
# Use the GP.chl dataset from previous testing block
test_that("merge_taxa() properly handles standard-cases", {
	expect_that(n1 <- merge_taxa(GP.chl, c("24341", "579085")), is_a("phyloseq"))
	expect_that(ntaxa(n1), equals(20L))
	# The first name is kept, others removed
	expect_that("579085" %in% taxa_names(n1), equals(FALSE))
	expect_that("24341" %in% taxa_names(n1),  equals(TRUE))
	# Try a 3-element merge
	expect_that(n2 <- merge_taxa(GP.chl, c("579085", "24341", "547579")), is_a("phyloseq"))
	expect_that(ntaxa(n2), equals(19L))
	# The first name is kept, others removed
	expect_that("579085" %in% taxa_names(n2), equals(TRUE))
	expect_that("24341"  %in% taxa_names(n2), equals(FALSE))
	expect_that("547579" %in% taxa_names(n2), equals(FALSE))	
	# Try again, but specify the retained OTU name as the 3rd one
	expect_that(n3 <- merge_taxa(GP.chl, c("579085", "24341", "547579"), "547579"), is_a("phyloseq"))
	# "547579" is kept, others removed
	expect_that("579085" %in% taxa_names(n3), equals(FALSE))
	expect_that("24341"  %in% taxa_names(n3), equals(FALSE))
	expect_that("547579" %in% taxa_names(n3), equals(TRUE))
	# Check that the remaining OTU has the sum of the values merged
	expect_that(get_sample(n3, "547579"), is_identical_to(
			get_sample(GP.chl, "24341") + get_sample(GP.chl, "547579") + get_sample(GP.chl, "579085")
		)
	)
})
test_that("merge_taxa() replaces disagreements in taxonomy with NA", {
	# Try a more difficult merge from a different subset
	GP20 <- prune_taxa(taxa_names(GlobalPatterns)[1:20], GlobalPatterns)
	
	# Arbitrary merge into taxa "951", NA in ranks after Phylum
	merge_these <- c("951", "586076", "141782", "30678", "30405")
	n5 <- merge_taxa(GP20, merge_these)
	# Test that none of the non-archetype taxa are left after the merge
	expect_that(all( !c("586076", "141782", "30678", "30405") %in% taxa_names(n5)), equals(TRUE))
	# Test that the archetype taxa remains
	expect_that( "951" %in% taxa_names(n5), equals(TRUE))
	# Test that the taxonomy is NA after Phylum
	n5_merged_taxonomy <- as(tax_table(n5), "matrix")["951", ]
	expect_that(any(is.na(n5_merged_taxonomy[1:2])), equals(FALSE))
	expect_that(all(is.na(n5_merged_taxonomy[3:7])), equals(TRUE))	
	
	# Test how well it works at a different level (say first or last ranks)
	merge_these <- c("1126", "31759")
	n6 <- merge_taxa(GP20, merge_these)
	# Test that the non-archetype taxa is gone
	expect_that( !"31759" %in% taxa_names(n6), equals(TRUE))
	# Test that the archetype taxa remains
	expect_that( "1126" %in% taxa_names(n6), equals(TRUE))
	# Test that the taxonomy is NA after Order
	n6_merged_taxonomy <- as(tax_table(n6), "matrix")["1126", ]
	expect_that( any(is.na(n6_merged_taxonomy[1:4])), equals(FALSE))
	expect_that( all(is.na(n6_merged_taxonomy[5:7])), equals(TRUE))	

	# Test that it works for differences at the first rank
	GP20f <- GP20
	tax_table(GP20f)[1, 1] <- "Bacteria"
	n7 <- merge_taxa(GP20f, taxa_names(GP20f)[1:2])
	# Should be all NA taxonomy
	expect_that( all(is.na(as(tax_table(n7), "matrix")[1, ])), equals(TRUE))

	# Test that it works for differences at the last rank
	# First, make the first taxa the same as "951"
	tax_table(GP20f)[1, ] <- tax_table(GP20f)["951", ]
	# Now change the last rank of this entry to something else
	tax_table(GP20f)[1, length(rank_names(GP20f))] <- "species_phyloseq_test"
	n8 <- merge_taxa(GP20f, c("951", taxa_names(GP20f)[1]))
	t951 <- as(tax_table(n8), "matrix")["951", ]	
	expect_that( sum(is.na(t951)), equals(1L))
	expect_that( is.na(t951[length(rank_names(n8))]), is_equivalent_to(TRUE))
	expect_that( t951[-7], is_identical_to(as(tax_table(GP20f), "matrix")["951", ][-7]))

	# Test that it works if the taxonomies completely agree
	GP20f <- GP20	
	# Make the first taxa the same as "951"
	tax_table(GP20f)[1, ] <- tax_table(GP20f)["951", ]
		
	merge_these <- c("549322", "951")
	n9   <- merge_taxa(GP20f, merge_these)
	n9t1 <- as(tax_table(n9), "matrix")["549322", ]
	# None should be NA
	expect_that(any(is.na(n9t1)), equals(FALSE))
	expect_that(length(n9t1), equals(7L))
	# Merge worked, "951" is gone.
	expect_that("951" %in% taxa_names(n9), equals(FALSE))	

})






