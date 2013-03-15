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
	expect_that(setequal(names(samobs), taxa_names(phy_tree(GP))[1:50]), is_false())
	expect_that(setequal(names(samobs), taxa_names(GP)[1:50]), is_false())
	expect_that(identical(names(samobs), taxa_names(GP)[1:50]), is_false())
})

# prune to just samobs OTUs
pGP = prune_taxa(names(samobs), GP)

test_that("The set of names should be the same after pruning, names(samobs)", {
	expect_that(setequal(names(samobs), taxa_names(phy_tree(pGP))), is_true())	
	expect_that(setequal(names(samobs), taxa_names(otu_table(pGP))), is_true())	
	expect_that(setequal(names(samobs), taxa_names(tax_table(pGP))), is_true())	
})

test_that("The set/order of taxa names after pruning should be consistent", {
	# set equal
	expect_that(setequal(taxa_names(pGP), taxa_names(phy_tree(pGP))), is_true())
	expect_that(setequal(taxa_names(otu_table(pGP)), taxa_names(phy_tree(pGP))), is_true())
	expect_that(setequal(taxa_names(tax_table(pGP)), taxa_names(phy_tree(pGP))), is_true())
	# identical
	expect_that(identical(taxa_names(pGP), taxa_names(phy_tree(pGP))), is_true())
	expect_that(identical(taxa_names(otu_table(pGP)), taxa_names(phy_tree(pGP))), is_true())
	expect_that(identical(taxa_names(tax_table(pGP)), taxa_names(phy_tree(pGP))), is_true())
	expect_that(identical(names(samobs), taxa_names(phy_tree(pGP))), is_false())
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
	expect_that(ent.logi, is_a("logical"))		
	ent.trim <- filter_taxa(enterotype, flist, TRUE)
	expect_that(ent.trim, is_a("phyloseq"))
	expect_that(sum(ent.logi), equals(ntaxa(ent.trim)))
	expect_that(prune_taxa(ent.logi, enterotype), is_identical_to(ent.trim)) 
	expect_that(ntaxa(ent.trim), equals(416L))
	expect_that(nsamples(ent.trim), equals(nsamples(enterotype)))	
})
