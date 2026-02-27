# load libraries
library("phyloseq"); library("testthat")
# # # # TESTS!
set.seed(888)

# Load GP dataset
data("GlobalPatterns")
GP <- GlobalPatterns
keepNames <- sample_names(GP)[5:7]

################################################################################
# prune_
################################################################################

test_that("Classes of pruned phyloseq and its components are as expected", {
	GP3 <- prune_samples(keepNames, GP)
	expect_identical(nsamples(GP3), 3L)
	expect_s4_class(GP3, "phyloseq")
	expect_s4_class(access(GP3, "sam_data"), "sample_data")
	expect_s4_class(access(GP3, "otu_table"), "otu_table")	
	expect_s3_class(access(GP3, "phy_tree"), "phylo")	
	expect_s4_class(access(GP3, "tax_table"), "taxonomyTable")		
	# Now try on instance without sample data (empty slot)
	GPnoSD  <- phyloseq(otu_table(GP), tax_table(GP))
	GP3noSD <- prune_samples(keepNames, GPnoSD)
	expect_identical(nsamples(GP3noSD), 3L)	
	expect_s4_class(access(GP3noSD, "otu_table"), "otu_table")	
	expect_null(access(GP3noSD, "sam_data"))
	expect_null(access(GP3noSD, "phy_tree"))	
	expect_s4_class(access(GP3noSD, "tax_table"), "taxonomyTable")	
})

test_that("prune_samples works on sample_data-only and otu_table-only data", {
	GPotu <- prune_samples(keepNames, access(GP, "otu_table", TRUE))
	GPsd  <- prune_samples(keepNames, access(GP, "sam_data", TRUE))	
	expect_identical(nsamples(GPotu), 3L)	
	expect_identical(nsamples(GPsd), 3L)			
	expect_s4_class(GPotu, "otu_table")
	expect_s4_class(GPsd, "sample_data")
	expect_identical(dim(GPotu), c(19216L, 3L))
	expect_identical(dim(GPsd), c(3L, 7L))			
})

# Coerce orientation for apply
if(taxa_are_rows(GP)){
  otumat = as(otu_table(GP), "matrix")
} else {
  otumat = t(as(otu_table(GP), "matrix"))
}
# Count in how many samples each OTU was observed more than 5 times.
samobs = apply(otumat, 1, function(x, m) sum(x > m), m=5L)
# Keep only the most prevalent 50 of these
samobs = sort(samobs, TRUE)[1:50]
# Shuffle the names on purpose.
samobs = sample(samobs, length(samobs), FALSE)

test_that("Initial order before pruning check is different", {
	expect_false(setequal(names(samobs), taxa_names(phy_tree(GP))[1:50]))
	expect_false(setequal(names(samobs), taxa_names(GP)[1:50]))
	expect_false(identical(names(samobs), taxa_names(GP)[1:50]))
})

# prune to just samobs OTUs
pGP = prune_taxa(names(samobs), GP)

test_that("The set of names should be the same after pruning, names(samobs)", {
	expect_true(setequal(names(samobs), taxa_names(phy_tree(pGP))))	
	expect_true(setequal(names(samobs), taxa_names(otu_table(pGP))))	
	expect_true(setequal(names(samobs), taxa_names(tax_table(pGP))))	
})

test_that("The set/order of taxa names after pruning should be consistent", {
	# set equal
	expect_true(setequal(taxa_names(pGP), taxa_names(phy_tree(pGP))))
	expect_true(setequal(taxa_names(otu_table(pGP)), taxa_names(phy_tree(pGP))))
	expect_true(setequal(taxa_names(tax_table(pGP)), taxa_names(phy_tree(pGP))))
	# identical
	expect_true(identical(taxa_names(pGP), taxa_names(phy_tree(pGP))))
	expect_true(identical(taxa_names(otu_table(pGP)), taxa_names(phy_tree(pGP))))
	expect_true(identical(taxa_names(tax_table(pGP)), taxa_names(phy_tree(pGP))))
	expect_false(identical(names(samobs), taxa_names(phy_tree(pGP))))
	# plot_tree(pGP, "sampledodge", nodeplotblank, label.tips="taxa_names", plot.margin=0.75)
})

## Add this as backup test
#' data("esophagus")
#' esophagus
#' plot(sort(taxa_sums(esophagus), TRUE), type="h", ylim=c(0, 50))
#' x1 = prune_taxa(taxa_sums(esophagus) > 10, esophagus) 
#' x2 = prune_taxa(names(sort(taxa_sums(esophagus), TRUE))[1:9], esophagus) 
#' identical(x1, x2)


################################################################################
# test filter_taxa and other filter methods.
################################################################################
library("genefilter")
data("enterotype")

test_that("filter_taxa gives correct, reliable logicals and pruning", {
	flist    <- filterfun(kOverA(5, 2e-05))
	ent.logi <- filter_taxa(enterotype, flist)
	expect_is(ent.logi, ("logical"))
	ent.trim <- filter_taxa(enterotype, flist, TRUE)
	expect_is(ent.trim, ("phyloseq"))
	expect_equal(sum(ent.logi), (ntaxa(ent.trim)))
	expect_identical(prune_taxa(ent.logi, enterotype), (ent.trim)) 
	expect_equal(ntaxa(ent.trim), (416L))
	expect_equal(nsamples(ent.trim), (nsamples(enterotype)))	
})