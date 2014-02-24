################################################################################
# Use testthat to test phyloseq transformation functions/methods
################################################################################
library("phyloseq"); library("testthat")
# # # # TESTS!
#set.seed(8888)

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
  expect_false(is.null(eso1), "eso1 is NULL, valid phyloseq construction failed.")
  expect_false(is.null(eso2), "eso2 is NULL, valid phyloseq construction failed.")
  expect_is(eso1, "phyloseq", "class of eso1 is not phyloseq")
  expect_is(eso2, "phyloseq", "class of eso2 is not phyloseq")
  expect_equal(ntaxa(esophagus), ntaxa(eso1), info="ntaxa eso1 doesn't match original after transformation.")
  expect_equal(ntaxa(esophagus), ntaxa(eso2), info="ntaxa eso2 doesn't match original after transformation.")
})

test_that("Test transform_sample_counts edge-cases", {
	data("esophagus")
	# Randomly pick the OTU and sample names to use for subsetting
	# (Will be different each time tests are run)
	OTUname1 = sample(taxa_names(esophagus), 1)
	samname1 = sample(sample_names(esophagus), 1)
	# Test that a one-OTU dataset still works
	# It throws a warning because the tree is being removed when it becomes just one tip.
	#suppressWarnings(eso1otu <- prune_taxa(sample(taxa_names(esophagus), 1), esophagus))
	expect_warning(try(eso1otu <- prune_taxa(OTUname1, esophagus), TRUE))
	try(eso1sam <- prune_samples(samname1, esophagus), TRUE)
	# Test eso1otu
	expect_equal(ntaxa(eso1otu), 1L)
	expect_equal(nsamples(eso1otu), 3L)
	# Behavior strange when tree removed by necessity.
	# In this case a "phyloseq" object with just an OTU table. 
	# Really, this should just be considered an object of class "otu_table".
	#expect_is(eso1otu, "phyloseq", "pruned-to-1-OTU esophagus not phyloseq")
	eso1otu = otu_table(eso1otu)
	expect_is(eso1otu, "otu_table")
	expect_equal(dim(otu_table(eso1otu)), c(1L, nsamples(esophagus)))	
	# Test eso1sam
	expect_is(eso1sam, "phyloseq", "pruned-to-1-sample esophagus not phyloseq")
	expect_is(otu_table(eso1sam), "otu_table",
						"pruned-to-1-sample esophagus OTU table not otu_table")
	expect_equal(nsamples(eso1sam), 1L)
	expect_equal(dim(otu_table(eso1sam)), c(ntaxa(esophagus), 1L))
	# Now test transform_sample_counts
	eso1samrank = eso1oturank = NULL
	try(eso1samrank  <- transform_sample_counts(eso1sam, rank), TRUE)
	try(eso1samrankt <- transform_sample_counts(t(eso1sam), rank), TRUE)
	try(eso1oturank  <- transform_sample_counts(eso1otu, rank), TRUE)
	try(eso1oturankt <- transform_sample_counts(t(eso1otu), rank), TRUE)
	expect_is(eso1samrank, "phyloseq")
	expect_is(eso1samrankt, "phyloseq")
	expect_is(eso1oturank, "otu_table")
	expect_is(eso1oturankt, "otu_table")
})
test_that("Test transform_sample_counts numerical result accuracy", {
	data("esophagus")
  es = esophagus
  # addition
  es1 = transform_sample_counts(es, function(x) x + 1)
  expect_equal(otu_table(es1), otu_table(es) + 1, tolerance=0.1, info="addition fail")
  # multiplication
	es1 = transform_sample_counts(es, function(x) x * 2.5)
	expect_equal(otu_table(es1), otu_table(es) * 2.5, tolerance=0.1, info="multiplication fail")
  # element-wise exponentiation
	es1 = transform_sample_counts(es, function(x) x ^ 2.5)
	expect_equal(otu_table(es1), otu_table(es) ^ 2.5, tolerance=0.1, info="exponentiation fail")
  # logarithm
	es1 = transform_sample_counts(es, function(x) log(x+10) )
	expect_equal(otu_table(es1), log(otu_table(es) + 10), tolerance=0.1, info="logarithm fail")  
  # Prune to a small subset. Need a test where "by-sample" matters. E.g. rank
  es  = prune_taxa(taxa_names(es)[1:5], esophagus)
	es1 = transform_sample_counts(es, rank)
  ans = c(5, 2, 4, 2, 2, 5, 2.5, 4, 2.5, 1, 5, 1.5, 1.5, 3.5, 3.5)
  ans = otu_table(matrix(ans, ntaxa(es), nsamples(es), FALSE,
                         list(taxa_names(es), sample_names(es))), taxa_are_rows=TRUE)
	expect_equal(otu_table(es1), ans, tolerance=0.1, info="rank fail")
  # test where "by-sample" matters, after transpose
	es  = prune_taxa(taxa_names(esophagus)[1:5], t(esophagus))
	es1 = transform_sample_counts(es, rank)
	ans = c(5, 2, 4, 2, 2, 5, 2.5, 4, 2.5, 1, 5, 1.5, 1.5, 3.5, 3.5)
	ans = otu_table(matrix(ans, nsamples(es), ntaxa(es), TRUE,
	                       list(sample_names(es), taxa_names(es))), taxa_are_rows=FALSE)
	expect_equal(otu_table(es1), ans, tolerance=0.1, info="rank fail")
})

