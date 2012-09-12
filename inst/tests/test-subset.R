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
	expect_that(access(GP3, "sam_data"),  is_a("sample_data"))
	expect_that(access(GP3, "otu_table"), is_a("otu_table"))	
	expect_that(access(GP3, "phy_tree"),      is_a("phylo"))	
	expect_that(access(GP3, "tax_table"),   is_a("taxonomyTable"))		
	# Now try on instance without sample data (empty slot)
	GPnoSD  <- phyloseq(otu_table(GP), tax_table(GP))
	GP3noSD <- prune_samples(keepNames, GPnoSD)
	expect_that(nsamples(GP3noSD), is_identical_to(3L))	
	expect_that(access(GP3noSD, "otu_table"), is_a("otu_table"))	
	expect_that(access(GP3noSD, "sam_data"),  is_a("NULL"))
	expect_that(access(GP3noSD, "phy_tree"),      is_a("NULL"))	
	expect_that(access(GP3noSD, "tax_table"),   is_a("taxonomyTable"))	
})

test_that("prune_samples works on sample_data-only and otu_table-only data", {
	GPotu <- prune_samples(keepNames, access(GP, "otu_table", TRUE))
	GPsd  <- prune_samples(keepNames, access(GP, "sam_data", TRUE))	
	expect_that(nsamples(GPotu), is_identical_to(3L))	
	expect_that(nsamples(GPsd), is_identical_to(3L))			
	expect_that(GPotu, is_a("otu_table"))
	expect_that(GPsd, is_a("sample_data"))
	expect_that(dim(GPotu), is_identical_to(c(19216L, 3L)))
	expect_that(dim(GPsd), is_identical_to(c(3L, 7L)))			
})

################################################################################
# test filter_taxa and other filter methods.
################################################################################
library("genefilter")
data("enterotype")

test_that("filter_taxa gives correct, reliable logicals and pruning", {
	flist    <- filterfun(kOverA(5, 2e-05))
	ent.logi <- filter_taxa(enterotype, flist)
	expect_that(ent.logi, is_a("logical"))		
	ent.trim <- filter_taxa(enterotype, flist, TRUE)
	expect_that(ent.trim, is_a("phyloseq"))
	expect_that(sum(ent.logi), equals(ntaxa(ent.trim)))
	expect_that(prune_taxa(ent.logi, enterotype), is_identical_to(ent.trim)) 
	expect_that(ntaxa(ent.trim), equals(416L))
	expect_that(nsamples(ent.trim), equals(nsamples(enterotype)))	
})
