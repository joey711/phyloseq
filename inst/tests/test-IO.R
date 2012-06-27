################################################################################
# Use testthat to test file import and resulting class (and values)
################################################################################
library("phyloseq"); library("testthat")
# # # # TESTS!

mothlist  <- system.file("extdata", "esophagus.fn.list", package="phyloseq")
mothgroup <- system.file("extdata", "esophagus.good.groups", package="phyloseq")
mothtree  <- system.file("extdata", "esophagus.tree", package="phyloseq")
cutoff    <- "0.10"
esophman  <- import_mothur(mothlist, mothgroup, mothtree, cutoff)	

test_that("import of esophagus dataset from mothur files in extdata/ produces a phyloseq object", {
	expect_that(esophman, is_a("phyloseq"))
})

test_that("The two phyloseq objects, example and just-imported, are identical", {
	data("esophagus")
	expect_that(esophagus, is_identical_to(esophman))
})

test_that("Test mothur file import on the (esophagus data).", {
	smlc <- show_mothur_list_cutoffs(mothlist)
	expect_that(smlc, is_identical_to(c("unique", "0.00", "0.01", "0.02", "0.03", "0.04", "0.05", "0.06", "0.07", "0.08", "0.09", "0.10")))	
})

test_that("abundances can be manipulated mathematically", {
	x1 <- as(otuTable(esophman), "matrix")
	expect_that(2*x1-x1, is_identical_to(x1) )
})

test_that("empty stuff is NULL", {
	expect_that(taxTab(esophman, FALSE), is_a("NULL"))
	expect_that(sampleData(esophman, FALSE), is_a("NULL"))
})

test_that("Expected classes of non-empty components", {
	expect_that(otuTable(esophman), is_a("otuTable"))
	expect_that(tre(esophman), is_a("phylo"))
})

test_that("imported files become S4 object", {
	expect_that(isS4(esophman), is_true())
})

test_that("show method output tests",{
	expect_that(esophman, prints_text("phyloseq-class experiment-level object"))
})


