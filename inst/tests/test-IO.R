################################################################################
# Use testthat to test file import and resulting class (and values)
################################################################################
library("phyloseq"); library("testthat")
# # # # TESTS!

################################################################################
# import_mothur tests
mothlist  <- system.file("extdata", "esophagus.fn.list.gz", package="phyloseq")
mothgroup <- system.file("extdata", "esophagus.good.groups.gz", package="phyloseq")
mothtree  <- system.file("extdata", "esophagus.tree.gz", package="phyloseq")
cutoff    <- "0.10"
esophman  <- import_mothur(mothlist, mothgroup, mothtree, cutoff)	

test_that("import_mothur: import of esophagus dataset from mothur files in extdata/ produces a phyloseq object", {
	expect_that(esophman, is_a("phyloseq"))
})

test_that("import_mothur: The two phyloseq objects, example and just-imported, are identical", {
	data("esophagus")
	expect_that(esophagus, is_identical_to(esophman))
})

test_that("import_mothur: Test mothur file import on the (esophagus data).", {
	smlc <- show_mothur_list_cutoffs(mothlist)
	expect_that(smlc, is_identical_to(c("unique", "0.00", "0.01", "0.02", "0.03", "0.04", "0.05", "0.06", "0.07", "0.08", "0.09", "0.10")))	
})

test_that("import_mothur: abundances can be manipulated mathematically", {
	x1 <- as(otuTable(esophman), "matrix")
	expect_that(2*x1-x1, is_identical_to(x1) )
})

test_that("import_mothur: empty stuff is NULL", {
	expect_that(taxTab(esophman, FALSE), is_a("NULL"))
	expect_that(sampleData(esophman, FALSE), is_a("NULL"))
})

test_that("import_mothur: Expected classes of non-empty components", {
	expect_that(otuTable(esophman), is_a("otuTable"))
	expect_that(tre(esophman), is_a("phylo"))
})

test_that("import_mothur: imported files become S4 object", {
	expect_that(isS4(esophman), is_true())
})

test_that("import_mothur: show method output tests",{
	expect_that(esophman, prints_text("phyloseq-class experiment-level object"))
})

################################################################################
# import_RDP tests
test_that("the import_RDP_otu function can properly read gzipped-example", {
	otufile <- system.file("extdata", "rformat_dist_0.03.txt.gz", package="phyloseq")
	ex_otu  <- import_RDP_otu(otufile)	

	expect_that(head(t(ex_otu)), prints_text("OTU Table:"))
	expect_that(ex_otu, is_a("otuTable"))
	expect_that(nspecies(ex_otu), equals(5276))
	expect_that(nsamples(ex_otu), equals(14))
	expect_that(sampleSums(ex_otu), is_a("numeric"))
})


################################################################################
# import_qiime tests
otufile <- system.file("extdata", "GP_otu_table_rand_short.txt.gz", package="phyloseq")
mapfile <- system.file("extdata", "master_map.txt", package="phyloseq")
trefile <- system.file("extdata", "GP_tree_rand_short.newick.gz", package="phyloseq")

t0 <- import_qiime(otufile, mapfile, trefile, showProgress=FALSE)
test_that("Class of import result is phyloseq-class", {
	expect_that(t0, is_a("phyloseq"))
})

test_that("Classes of components are as expected", {
	expect_that(otuTable(t0), is_a("otuTable"))
	expect_that(taxTab(t0), is_a("taxonomyTable"))
	expect_that(samData(t0), is_a("sampleData"))
	expect_that(tre(t0), is_a("phylo"))		
})

test_that("Changing the chunk.size does not affect resulting tables", {
	t1 <- import_qiime(otufile, mapfile, trefile, chunk.size=300L, showProgress=FALSE)
	t2 <- import_qiime(otufile, mapfile, trefile, chunk.size=13L, showProgress=FALSE)
	expect_that(t0, is_identical_to(t1))
	expect_that(t1, is_identical_to(t2))
})	

test_that("Features of the abundance data are consistent, match known values", {
	expect_that(sum(speciesSums(t0)), equals(1269671L))
	expect_that(sum(speciesSums(t0)==0), equals(5L))
	expect_that(sum(speciesSums(t0)>=100), equals(183L))
	expect_that(sum(speciesSums(t0)), equals(sum(sampleSums(t0))))
	expect_that(sum(sampleSums(t0) > 10000L), equals(20L))
	expect_that(nsamples(t0), equals(26L))
	expect_that(nspecies(t0), equals(500L))
})

################################################################################