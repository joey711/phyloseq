################################################################################
# Use testthat to test that distance methods return correct results.
################################################################################
library("phyloseq")
library("testthat")
################################################################################
# UniFrac testing section. Relies on pre-computed results from pycogent
# The relevant python code is saved in `extdata/gp500-pycogent.py`
################################################################################
# Define the test-example phyloseq dataset. A random subsample from Global Patterns
treeFile = system.file("extdata", "GP_tree_rand_short.newick.gz", package="phyloseq")
GP500File = system.file("extdata", "GP_otu_table_rand_short.txt.gz", package = "phyloseq")
GP500 = import_qiime(GP500File, treefilename = treeFile)
# # Example if you want to re-create the test files for calculating with pyCogent
#export_env_file(GP500, file = "~/Downloads/gp500test.env.txt", writeTree = FALSE)
#ape::write.tree(phy_tree(GP500), file = "~/Downloads/gp500test.tree")
# Now import the results with read.table()
gp500_uuf = read.csv(system.file("extdata", "gp500-uuf.csv", package = "phyloseq"), header = FALSE, fill = TRUE)
gp500_wuf = read.csv(system.file("extdata", "gp500-wuf.csv", package = "phyloseq"), header = FALSE, fill = TRUE)
gp500_wufu = read.csv(system.file("extdata", "gp500-wufu.csv", package = "phyloseq"), header = FALSE, fill = TRUE)
# Add the sample names
colnames(gp500_uuf) <- rownames(gp500_uuf) <- colnames(gp500_wuf) <- rownames(gp500_wuf) <- colnames(gp500_wufu) <- rownames(gp500_wufu) <- sample_names(GP500)
# Coerce to Distance Matrices for comparison `"dist"` class
gp500_wufu <- as.dist(gp500_wufu)
gp500_wuf <- as.dist(gp500_wuf)
gp500_uuf <- as.dist(gp500_uuf)
# Define numerical tolerance
tol = 0.00000001
test_that("UniFrac produces correct values on an example subset from Global Patterns. 'Correct' values are results from pyCogent", {
  # Using UniFrac function directly
  expect_equal(gp500_wufu, UniFrac(GP500, weighted = TRUE, normalized = FALSE), check.attributes = FALSE, tolerance = tol,
               label = "`UniFrac`: Weighted but Un-normalized UniFrac results did not match reference answer.")
  expect_equal(gp500_wuf, UniFrac(GP500, weighted = TRUE), check.attributes = FALSE, tolerance = tol,
               label = "`UniFrac`: Weighted, normalized UniFrac results did not match reference answer.")
  expect_equal(gp500_uuf, UniFrac(GP500, weighted = FALSE), check.attributes = FALSE, tolerance = tol,
               label = "`UniFrac`: Unweighted UniFrac results did not match reference answer.")
  # Using the `distance` wrapper
  expect_equal(gp500_wufu, distance(GP500, "unifrac", weighted = TRUE, normalized = FALSE), check.attributes = FALSE, tolerance = tol,
               label = "`distance`: Weighted but Un-normalized UniFrac results did not match reference answer.")
  expect_equal(gp500_wuf, distance(GP500, "unifrac", weighted = TRUE), check.attributes = FALSE, tolerance = tol,
               label = "`distance`: Weighted, normalized UniFrac results did not match reference answer.")
  expect_equal(gp500_uuf, distance(GP500, "unifrac", weighted = FALSE), check.attributes = FALSE, tolerance = tol,
               label = "`distance`: Unweighted UniFrac results did not match reference answer.")  
  # Make sure reference files are different (at the very least)
  expect_false({isTRUE(all.equal(gp500_uuf, gp500_wuf, check.attributes = FALSE, tolerance = 0.01))},
               label = "The reference matrices for UniFrac testing should be different, but were not. uuf/wuf")
  expect_false({isTRUE(all.equal(gp500_uuf, gp500_wufu, check.attributes = FALSE, tolerance = 0.01))},
               label = "The reference matrices for UniFrac testing should be different, but were not. uuf/wufu")
  expect_identical(distance(GP500, "wunifrac"), distance(GP500, "unifrac", weighted = TRUE),
                   label = "wunifrac output is not identical to unifrac with weighted=T flag")
})
test_that("Check that regular-expression matching for unifrac method flag is working", {
  expect_identical(distance(GP500, "w-UniFrac"), distance(GP500, "unifrac", weighted = TRUE))
  expect_identical(distance(GP500, "weighted-UniFrac"), distance(GP500, "unifrac", weighted = TRUE))
  expect_identical(distance(GP500, "unweighted-UniFrac"), distance(GP500, "unifrac"))
  expect_identical(distance(GP500, "u-UniFrac"), distance(GP500, "unifrac"))
})
################################################################################
# Test other distances against their expected dispatch explicit calculation
################################################################################
test_that("Test accurate dispatch for other distances", {
  otumatgp500 = t(as(otu_table(GP500), "matrix"))
  # Test for all vegdist, phyloseq object
  expect_equal(
    lapply(distanceMethodList$vegdist, vegan::vegdist, x = otumatgp500),
    lapply(distanceMethodList$vegdist, distance, physeq = GP500),
    check.attributes = FALSE, tolerance = tol)
  # Test for all vegdist, OTU table
  expect_equal(
    lapply(distanceMethodList$vegdist, vegan::vegdist, x = otumatgp500),
    lapply(distanceMethodList$vegdist, distance, physeq = otu_table(GP500)),
    check.attributes = FALSE, tolerance = tol)
  # Test for all betadiver, phyloseq object
  expect_equal(
    lapply(distanceMethodList$betadiver, vegan::betadiver, x = otumatgp500),
    lapply(distanceMethodList$betadiver, distance, physeq = GP500),
    check.attributes = FALSE, tolerance = tol)
  # Test for all betadiver, OTU table
  expect_equal(
    lapply(distanceMethodList$betadiver, vegan::betadiver, x = otumatgp500),
    lapply(distanceMethodList$betadiver, distance, physeq = otu_table(GP500)),
    check.attributes = FALSE, tolerance = tol) 
  # Test for all dist, phyloseq object
  expect_equal(
    lapply(distanceMethodList$dist, stats::dist, x = otumatgp500),
    lapply(distanceMethodList$dist, distance, physeq = GP500),
    check.attributes = FALSE, tolerance = tol)
  # Test for all dist, OTU table
  expect_equal(
    lapply(distanceMethodList$dist, stats::dist, x = otumatgp500),
    lapply(distanceMethodList$dist, distance, physeq = otu_table(GP500)),
    check.attributes = FALSE, tolerance = tol) 
  
  # DPCoA
  #"dpcoa"
  # phyloseq object
  # TOO SLOW to test routinely. Commenting out.
  #   expect_equal(
  #     as.dist(DPCoA(GP500)$RaoDis, diag=FALSE),
  #     distance(GP500, "dpcoa"),
  #     check.attributes = FALSE, tolerance = tol)
  # OTU table doesn't have a tree
  expect_error(DPCoA(otu_table(GP500)))
  
  # JSD
  #"jsd"
  # phyloseq object
  expect_equal(
    phyloseq:::JSD(GP500),
    distance(GP500, "jsd"),
    check.attributes = FALSE, tolerance = tol)
  # OTU table
  expect_equal(
    phyloseq:::JSD(otu_table(GP500)),
    distance(otu_table(GP500), "jsd"),
    check.attributes = FALSE, tolerance = tol)
})
################################################################################
