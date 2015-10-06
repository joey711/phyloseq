# testthat tests don't do anything when successful.
library("phyloseq"); packageVersion("phyloseq")
library("testthat"); packageVersion("testthat")


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
test_that("merge_phyloseq: Break apart GP based on human-association,
          then merge back together.", {
	data(GlobalPatterns)
	GP  <- prune_taxa(taxa_names(GlobalPatterns)[1:100], GlobalPatterns)
	sample_data(GP)$human <- factor(get_variable(GP, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue"))
	h1 <- subset_samples(GP, human=="TRUE")
	h0 <- subset_samples(GP, human=="FALSE")
	GP1 <- merge_phyloseq(h0, h1)

	# The species order is fixed to the tree, 
	# so should be the same between the original and merged
	expect_that(taxa_names(GP), is_identical_to(taxa_names(GP1)))
	expect_that(phy_tree(h1), is_identical_to(phy_tree(h0)))

	# However, the sample order has been shuffled by the split/merge. 
	# Fix the sample order by re-ordering the otu_table, and reassigning
	sa.order <- sample_names(GP)
	sa.order <- sa.order[sa.order %in% sample_names(GP1)]
	otu_table(GP1) <- otu_table(GP1)[, sa.order]
	expect_equal(sample_names(GP), sample_names(GP1))
	expect_equal(sample_names(sample_data(GP)), sample_names(sample_data(GP1)))
	# Sample data entries are the same, irrespective of factor levels
	GPfactors = which(sapply(sample_data(GP1), inherits, "factor"))
	for(j in GPfactors){
	  expect_equal(as.character(get_variable(GP, j)),
	               as.character(get_variable(GP1, j)))
	}
	# Reconcile factor level order for remaining tests
	GP1factors = which(sapply(sample_data(GP1), inherits, "factor"))
	for(j in names(GP1factors)){
	  varj = as.character(get_variable(GP1, j))
	  sample_data(GP1)[, j] <- factor(varj, levels = sort(unique(varj)))
	}
	GPfactors = which(sapply(sample_data(GP), inherits, "factor"))
	for(j in names(GPfactors)){
	  varj = as.character(get_variable(GP, j))
	  sample_data(GP)[, j] <- factor(varj, levels = sort(unique(varj)))
	}
	# Check a specific variable
	expect_equal(sample_data(GP1)$SampleType,
	             sample_data(GP)$SampleType)
	# Should be fixed now. Full object and components now identical
	expect_equal(GP1, GP) 
	expect_identical(GP1, GP)
	expect_that(otu_table(GP1), is_identical_to(otu_table(GP)))
	expect_that(tax_table(GP1), is_identical_to(tax_table(GP)))
	expect_that(phy_tree(GP1), is_identical_to(phy_tree(GP)))
	
	## Check factor levels
	# The set
	expect_identical(sort(levels(sample_data(GP1)$SampleType)),
	                 sort(levels(sample_data(GP)$SampleType)))
	# The order
	expect_identical(levels(sample_data(GP1)$SampleType),
	                 levels(sample_data(GP)$SampleType))	
	# Overall
	expect_identical(sample_data(GP1), sample_data(GP))
	expect_identical(droplevels(sample_data(GP1)), droplevels(sample_data(GP)))
	
	# Check variable names are all there (set)
	expect_equal(
	  object = sort(intersect(colnames(sample_data(GP1)), colnames(sample_data(GP)))),
	  expected = sort(colnames(sample_data(GP1))))
	# Check column classes
	expect_equal(sapply(sample_data(GP1), class), sapply(sample_data(GP), class))
	# Check column names
	expect_equal(colnames(sample_data(GP1)), colnames(sample_data(GP)))
	# Check sample name order
	expect_equal(sample_names(sample_data(GP1)), sample_names(sample_data(GP)))
	expect_equal(sample_names(GP1), sample_names(GP))
})

################################################################################
# tax_glom
# Load data
data("GlobalPatterns")
GP.chl = subset_taxa(GlobalPatterns, Phylum == "Chlamydiae")
test_that("the tax_table slot is identical whether tax_glom()ed by itself or as component", {
	expect_that(tax_glom(tax_table(GP.chl), "Family"), is_a("taxonomyTable"))
	expect_that(n1<-tax_glom(GP.chl, "Family"), is_a("phyloseq"))
	expect_that(ntaxa(n1), equals(4L))
	expect_that(
		tax_glom(tax_table(GP.chl), taxrank="Family"),
		is_equivalent_to(tax_table(tax_glom(GP.chl, taxrank="Family")))
	)
	n1 = as(tax_glom(tax_table(GP.chl), taxrank="Family", NArm=FALSE), "matrix")[, "Family"]
	n2 = tax_glom(GP.chl, taxrank="Family", NArm=FALSE)
  expect_true(setequal(n1, as(tax_table(n2), "matrix")[, "Family"]))
	expect_that(ntaxa(n2), equals(5L))	
})
test_that("tax_glom() handles clearly agglomeration to one taxa", {
	expect_that(n1 <- tax_glom(GP.chl, "Phylum"), gives_warning())
	expect_that(n1, is_a("phyloseq"))
	expect_that(ntaxa(n1), equals(1L))
	expect_that(access(n1, "phy_tree"), is_a("NULL"))
})
test_that("tax_glom() can handle even the highest rank glom", {
  expect_warning(tax_glom(GP.chl, "Kingdom"))
  gpk = tax_glom(GlobalPatterns, "Kingdom")
  expect_is(gpk, "phyloseq")
  expect_equivalent(ntaxa(gpk), 2)
  expect_equivalent(taxa_sums(gpk), c(195598, 28021080))
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
	# Try a 3-element merge, check that the largest-count remains.
  OTUIDs = c("579085", "24341", "547579")
	biggestOTU = names(which.max(taxa_sums(GP.chl)[OTUIDs]))
  # Perform the merge of `OTUIDs`, and check the resulting class while at it.
	expect_is(n2 <- merge_taxa(GP.chl, OTUIDs), "phyloseq")
  # Check that there are now the correct, fewer number of OTUs
	expect_equal(ntaxa(n2), (ntaxa(GP.chl)-length(OTUIDs)+1))
	# The biggest OTU is kept, others merged
  expect_true(biggestOTU %in% taxa_names(n2))
	expect_true(!any(setdiff(OTUIDs, biggestOTU) %in% taxa_names(n2)))
  # Merge again, but only use the tax_table. No counts changes default retained to first in vector
	expect_is(n2b <- merge_taxa(tax_table(GP.chl), OTUIDs), "taxonomyTable")
	# Check that there are now the correct, fewer number of OTUs
	expect_equal(ntaxa(n2b), (ntaxa(GP.chl)-length(OTUIDs)+1))
	# The biggest OTU is kept, others merged
	expect_true(OTUIDs[1] %in% taxa_names(n2b))
	expect_true(!any(setdiff(OTUIDs, OTUIDs[1]) %in% taxa_names(n2b)))
	# Merge again, but specify the retained OTU name as the 3rd one, rather than the default
	expect_that(n3 <- merge_taxa(GP.chl, eqtaxa=OTUIDs, archetype=OTUIDs[3]), is_a("phyloseq"))
	# "547579" is kept, others removed
	expect_true(OTUIDs[3] %in% taxa_names(n3))
	expect_true(!any(setdiff(OTUIDs, OTUIDs[3]) %in% taxa_names(n3)))
	# Check that the remaining OTU has the sum of the values merged
	expect_identical(get_sample(n3, OTUIDs[3]), 
    colSums(as(otu_table(GP.chl), "matrix")[OTUIDs, ]))
})
test_that("merge_taxa() replaces disagreements in taxonomy with NA", {
	# Try a more difficult merge from a different subset
	GP20 <- prune_taxa(taxa_names(GlobalPatterns)[1:20], GlobalPatterns)
	# Arbitrary merge into taxa "951", NA in ranks after Phylum
	OTUIDs = c("951", "586076", "141782", "30678", "30405")
	biggestOTU = names(which.max(taxa_sums(GP20)[OTUIDs]))  
	n5 = merge_taxa(GP20, OTUIDs)
	# The biggest OTU is kept, others merged
	expect_true(biggestOTU %in% taxa_names(n5))
	expect_true(!any(setdiff(OTUIDs, biggestOTU) %in% taxa_names(n5)))
	# The taxonomy should be NA_character_ after Phylum (OTUIDs chosen carefully in this case)
	n5_merged_taxonomy <- as(tax_table(n5), "matrix")[biggestOTU, ]
	expect_true(!any(is.na(n5_merged_taxonomy[1:2])))
	expect_true(all(is.na(n5_merged_taxonomy[3:7])))	
	# Test how well it works at a different level (say first or last ranks)
	OTUIDs <- c("1126", "31759")
	biggestOTU = names(which.max(taxa_sums(GP20)[OTUIDs]))  
	n6 <- merge_taxa(GP20, OTUIDs)
	# The biggest OTU is kept, others merged
	expect_true(biggestOTU %in% taxa_names(n6))
	expect_true(!any(setdiff(OTUIDs, biggestOTU) %in% taxa_names(n6)))
	# Test that the taxonomy is NA after Order
	n6_merged_taxonomy <- as(tax_table(n6), "matrix")[biggestOTU, ]
	expect_true( !any(is.na(n6_merged_taxonomy[1:4])) )
	expect_true( all(is.na(n6_merged_taxonomy[5:7])) )	
	# Test that it works for differences at the first rank
	GP20f <- GP20
	tax_table(GP20f)[1, 1] <- "Bacteria"
  OTUIDs = taxa_names(GP20f)[1:2]
	biggestOTU = names(which.max(taxa_sums(GP20f)[OTUIDs]))  
  expect_is(n7 <- merge_taxa(GP20f, OTUIDs), "phyloseq")
	# Should be all NA taxonomy
	expect_that( all(is.na(as(tax_table(n7), "matrix")[biggestOTU, ])), equals(TRUE))
	# Test that it works for differences at the last rank
	# First, make the first taxa the same as "951"
	tax_table(GP20f)[1, ] <- tax_table(GP20f)["951", ]
	# Now change the last rank of this entry to something else
	tax_table(GP20f)[1, length(rank_names(GP20f))] <- "species_phyloseq_test"
  OTUIDs = c("951", biggestOTU)
	biggestOTU = names(which.max(taxa_sums(GP20f)[OTUIDs]))
	expect_is(n8 <- merge_taxa(GP20f, OTUIDs), "phyloseq")
	t951 <- as(tax_table(n8), "matrix")[biggestOTU, ]	
	expect_equal( sum(is.na(t951)), 1L )
	expect_true( is.na(t951[length(rank_names(n8))]) )
	expect_identical( t951[-7],  as(tax_table(GP20f), "matrix")["951", ][-7] )
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
test_that("merge_taxa() properly handles different types and orders of taxa specified by the eqtaxa and archetype arguments, and also handles refseq data", {
  # Test merge_taxa on data with a reference sequence file.
  otufile <- system.file("extdata", "GP_otu_table_rand_short.txt.gz", package="phyloseq")
  mapfile <- system.file("extdata", "master_map.txt", package="phyloseq")
  trefile <- system.file("extdata", "GP_tree_rand_short.newick.gz", package="phyloseq")
  rs_file <- system.file("extdata", "qiime500-refseq.fasta", package="phyloseq")
  rs0 <- import_qiime(otufile, mapfile, trefile, rs_file)
  rs1 = merge_taxa(rs0, c("71074", "10517", "8096"))
  rs2 = merge_taxa(rs0, c("71074", "8096", "10517"), "71074")
  rs3 = merge_taxa(rs0, c("71074", "10517", "8096"), 3)
  rs4 = merge_taxa(rs0, c("8096", "71074", "10517"))
  # rs1 and rs2 should be identical
  # rs3 and rs4 should be identical
  expect_equivalent(rs1, rs2)
  expect_true(!identical(rs1, rs3))
  expect_equivalent(rs3, rs4)
  # double-check that components are all there
  expect_that(length(getslots.phyloseq(rs1)), equals(5L))
  expect_that(length(getslots.phyloseq(rs2)), equals(5L))
  expect_that(length(getslots.phyloseq(rs3)), equals(5L))
  expect_that(length(getslots.phyloseq(rs4)), equals(5L))
  # The number of taxa should be the same as the original less two
  expect_that(ntaxa(rs1), equals(ntaxa(rs0)-2L))
  expect_that(ntaxa(rs2), equals(ntaxa(rs0)-2L))
  expect_that(ntaxa(rs3), equals(ntaxa(rs0)-2L))
  expect_that(ntaxa(rs4), equals(ntaxa(rs0)-2L))	
  # merge_taxa() errors when a bad archetype is provided
  # Throws error because keepIndex is NULL
  expect_that(merge_taxa(rs0, c("71074", "10517", "8096"), "wtf"), throws_error())
  # Throws error because keepIndex is not part of eqtaxa (logic error, invalid merge)
  expect_that(merge_taxa(rs0, c("71074", "10517", "8096"), "13662"), throws_error())
})
