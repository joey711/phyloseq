transform_test_table <- function(n_taxa=3L, n_samples=2L, taxa_rows=TRUE) {
  counts <- matrix(rep(c(0, 2, 6, 4, 0, 8), length.out=n_taxa * n_samples),
                   nrow=n_taxa,
                   dimnames=list(paste0("taxon", seq_len(n_taxa)),
                                 paste0("sample", seq_len(n_samples))))
  if (!taxa_rows) counts <- t(counts)
  otu_table(counts, taxa_are_rows=taxa_rows)
}

test_that("sample transformations preserve orientation, names and singleton dimensions", {
  for (taxa_rows in c(TRUE, FALSE)) {
    for (n_taxa in c(1L, 3L)) {
      for (n_samples in c(1L, 2L)) {
        input <- transform_test_table(n_taxa, n_samples, taxa_rows)
        original <- input
        result <- transform_sample_counts(input, function(x) unname(log1p(x)))

        expect_s4_class(result, "otu_table")
        expect_identical(taxa_are_rows(result), taxa_rows)
        expect_identical(dim(result), dim(input))
        expect_identical(dimnames(result), dimnames(input))
        expect_equal(as(result, "matrix"), log1p(as(input, "matrix")))
        expect_identical(input, original)
      }
    }
  }
})

test_that("extra arguments reach the transformation function", {
  for (taxa_rows in c(TRUE, FALSE)) {
    input <- transform_test_table(taxa_rows=taxa_rows)
    result <- transform_sample_counts(input, function(x, offset) x + offset, offset=7)
    expect_equal(as(result, "matrix"), as(input, "matrix") + 7)
  }
})

test_that("callbacks run once per actual sample with the original taxon names", {
  for (taxa_rows in c(TRUE, FALSE)) {
    input <- transform_test_table(taxa_rows=taxa_rows)
    seen <- list()
    result <- transform_sample_counts(input, function(x) {
      seen[[length(seen) + 1L]] <<- x
      x / sum(x)
    })

    expect_length(seen, nsamples(input))
    counts <- as(if (taxa_rows) input else t(input), "matrix")
    for (i in seq_len(nsamples(input))) {
      expect_identical(seen[[i]], counts[, i])
    }
    expect_equal(unname(sample_sums(result)), rep(1, nsamples(input)))
  }
})

test_that("transformations preserve every other phyloseq component", {
  counts <- transform_test_table()
  samples <- sample_data(data.frame(group=c("control", "case"),
                                    row.names=sample_names(counts)))
  taxonomy <- tax_table(matrix(c("A", "B", "C"), ncol=1,
                               dimnames=list(taxa_names(counts), "Genus")))
  tree <- ape::read.tree(text="(taxon1:1,(taxon2:1,taxon3:1):1);")
  sequences <- Biostrings::DNAStringSet(setNames(c("AAA", "CCC", "GGG"),
                                               taxa_names(counts)))
  input <- phyloseq(counts, samples, taxonomy, tree, sequences)

  for (taxa_rows in c(TRUE, FALSE)) {
    current <- if (taxa_rows) input else t(input)
    original <- current
    result <- transform_sample_counts(current, log1p)

    expect_s4_class(result, "phyloseq")
    expect_identical(taxa_are_rows(result), taxa_rows)
    expect_equal(as(otu_table(result), "matrix"), log1p(as(otu_table(current), "matrix")))
    expect_identical(sample_data(result), sample_data(current))
    expect_identical(tax_table(result), tax_table(current))
    expect_identical(phy_tree(result), phy_tree(current))
    expect_identical(refseq(result), refseq(current))
    expect_identical(current, original)
  }
})

test_that("integer and double transformations are both supported", {
  input <- transform_test_table()
  for (fun in list(rank, function(x) as.integer(x > 0))) {
    result <- transform_sample_counts(input, fun)
    expected <- apply(as(input, "matrix"), 2, fun)
    dimnames(expected) <- dimnames(input)
    expect_identical(as(result, "matrix"), expected)
  }
})

test_that("invalid callback results fail before replacing the OTU table", {
  input <- transform_test_table()
  original <- input
  for (fun in list(sum, function(x) x[-1], function(x) as.character(x),
                   as.list, function(x) matrix(x, ncol=1),
                   function(x) if (x[1] == 0) x else x[-1])) {
    expect_error(transform_sample_counts(input, fun), "numeric vector.*same length")
    expect_identical(input, original)
  }
  expect_error(transform_sample_counts(input, 1), "`fun` must be a function")
})

test_that("the legacy transformation alias has the same behavior", {
  input <- transform_test_table()
  expect_identical(transformSampleCounts(input, log1p),
                   transform_sample_counts(input, log1p))
})

test_that("threshold rank closures honor zero handling and rank arguments", {
  x <- c(0, 5, 2, 0)
  expect_equal(threshrankfun(2, keep0s=TRUE, ties.method="min")(x), c(0, 4, 3, 0))
  for (keep_zeros in c(TRUE, FALSE)) {
    expect_identical(threshrankfun(2, keep0s=keep_zeros, ties.method="min")(x),
                     threshrank(x, 2, keep0s=keep_zeros, ties.method="min"))
  }
})
