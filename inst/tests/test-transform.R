################################################################################
# Use testthat to test phyloseq transformation functions/methods
################################################################################
library("phyloseq"); library("testthat")
# # # # TESTS!
set.seed(8888)

################################################################################
test_that("Can transform_sample_counts of an OTU table that is either orientation", {
  data("esophagus")
  OTU0 = otu_table(esophagus)
  OTU1 = transform_sample_counts(OTU0, rank)
  OTU2 = transform_sample_counts(t(OTU0), rank)
  expect_that(identical(ntaxa(OTU0), ntaxa(OTU1)), is_true(),
              "ntaxa OTU1 doesn't match original after transformation.")
  expect_that(identical(ntaxa(OTU0), ntaxa(OTU2)), is_true(),
              "ntaxa OTU2 doesn't match original after transformation.")
})
test_that("Can transform_sample_counts on phyloseq with either orientation", {
  data("esophagus")
  eso1 = eso2 = NULL
  try(eso1 <- transform_sample_counts(esophagus, rank), TRUE)
  try(eso2 <- transform_sample_counts(t(esophagus), rank), TRUE)
  expect_that(is.null(eso1), is_false(),
              "eso1 is NULL, valid phyloseq construction failed.")
  expect_that(is.null(eso2), is_false(),
              "eso2 is NULL, valid phyloseq construction failed.")
  expect_that(class(eso1), is_equivalent_to("phyloseq"),
              "class of eso1 is not phyloseq")
  expect_that(class(eso2), is_equivalent_to("phyloseq"),
              "class of eso2 is not phyloseq")
  expect_that(identical(ntaxa(esophagus), ntaxa(eso1)), is_true(),
              "ntaxa eso1 doesn't match original after transformation.")
  expect_that(identical(ntaxa(esophagus), ntaxa(eso2)), is_true(),
              "ntaxa eso2 doesn't match original after transformation.")
})