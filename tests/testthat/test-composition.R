composition_fixture <- function(zeros=FALSE, taxa_rows=TRUE) {
  values <- matrix(c(1, 2, 4, 8, 2, 2, 1, 3, 6, 3, 2, 1), nrow=4,
                    dimnames=list(paste0("t", 1:4), c("01", "s2", "s3")))
  if (zeros) values[1, 1] <- 0
  counts <- otu_table(if (taxa_rows) values else t(values), taxa_rows)
  taxonomy <- tax_table(matrix(c("A", "A", "B", "B"), ncol=1,
                               dimnames=list(rownames(values), "Genus")))
  samples <- sample_data(data.frame(group=factor(c("case", "control", "case")),
                                    age=c(30, 40, 50), row.names=colnames(values)))
  tree <- ape::read.tree(text="((t1:1,t2:1):1,(t3:1,t4:1):1);")
  sequences <- Biostrings::DNAStringSet(setNames(c("AAA", "CCC", "GGG", "TTT"),
                                                rownames(values)))
  phyloseq(counts, taxonomy, samples, tree, sequences)
}

test_that("composition checks distinguish scales and report invalid samples", {
  values <- matrix(c(1, 2, 4, 0, NA, -1), 3,
                    dimnames=list(c("a", "b", "c"), c("good", "bad")))
  report <- check_composition(values, scale="counts")
  expect_false(report$valid)
  expect_identical(report$samples$Valid, c(TRUE, FALSE))
  expect_equal(report$samples$NonFinite, c(0, 1))
  expect_equal(report$samples$Negative, c(0, 1))
  expect_error(check_composition(values, stop_on_error=TRUE), "Non-finite")
  for (bad in c(NA_real_, NaN, Inf, -Inf, -1)) {
    current <- values[, 1, drop=FALSE]
    current[1, 1] <- bad
    expect_false(check_composition(current)$valid)
  }
  values <- matrix(c(0.25, 0.75), 2, dimnames=list(c("a", "b"), "s"))
  expect_true(check_composition(values, scale="proportion")$valid)
  expect_false(check_composition(values, scale="counts")$valid)
  expect_true(check_composition(values, scale="abundance")$valid)
  expect_false(check_composition(values * 2, scale="proportion")$valid)
  values[,] <- c(-1, 1)
  expect_true(check_composition(values, scale="clr")$valid)
  expect_false(check_composition(values, scale="abundance")$valid)
  expect_false(check_composition(values + 0.01, scale="clr")$valid)
})

test_that("empty samples and ambiguous identifiers are never silently accepted", {
  values <- matrix(c(0, 1, 0, 0), 2,
                    dimnames=list(c("never", "present"), c("ok", "empty")))
  report <- check_composition(values)
  expect_false(report$valid)
  expect_identical(report$samples$Valid, c(TRUE, FALSE))
  expect_identical(report$zero_taxa, "never")
  expect_identical(report$samples$Sample, c("ok", "empty"))
  rownames(values) <- c("a", "a")
  expect_match(paste(check_composition(values)$issues, collapse=" "), "unique.*taxon")
  rownames(values) <- c("a", "b")
  colnames(values) <- c("", "s2")
  expect_match(paste(check_composition(values)$issues, collapse=" "), "unique.*sample")
  dimnames(values) <- NULL
  expect_false(check_composition(values)$valid)
  expect_error(check_composition(matrix(1+1i, 2, 2)), "real numeric")
  expect_error(check_composition(matrix(numeric(), 0, 1)), "non-empty")
  expect_error(check_composition(values, tolerance=0), "tolerance")
  expect_error(check_composition(values, stop_on_error=NA), "stop_on_error")
})

test_that("TreeSummarizedExperiment conversion preserves all components and orientations", {
  skip_if_not_installed("TreeSummarizedExperiment", "2.14.0")
  for (taxa_rows in c(TRUE, FALSE)) {
    original <- composition_fixture(taxa_rows=taxa_rows)
    tse <- as_tree_experiment(original)
    counts <- as(otu_table(original), "matrix")
    if (!taxa_rows) counts <- t(counts)
    expect_s4_class(tse, "TreeSummarizedExperiment")
    expect_identical(SummarizedExperiment::assay(tse, "counts"), counts)
    expect_identical(as.data.frame(SummarizedExperiment::colData(tse)),
                     as(sample_data(original), "data.frame"))
    expect_identical(as.matrix(SummarizedExperiment::rowData(tse)), as(tax_table(original), "matrix"))
    expect_identical(TreeSummarizedExperiment::rowTree(tse), phy_tree(original))
    expect_identical(TreeSummarizedExperiment::referenceSeq(tse), refseq(original))
    expect_identical(sample_names(original), colnames(tse))
    expect_true(methods::validObject(tse))
  }
  input <- otu_table(matrix(2, dimnames=list("taxon", "01")), TRUE)
  tse <- as_tree_experiment(input, assay="original")
  expect_identical(dim(tse), c(1L, 1L))
  expect_identical(dimnames(tse), dimnames(input))
  expect_null(TreeSummarizedExperiment::rowTree(tse))
  expect_error(as_tree_experiment(matrix(1, 1, 1)), "phyloseq or otu_table")
})

test_that("new assays retain raw data, metadata and existing transformation records", {
  skip_if_not_installed("TreeSummarizedExperiment", "2.14.0")
  input <- composition_fixture()
  original <- input
  tse <- transform_composition(input, "relative", input_scale="counts")
  tse <- transform_composition(tse, "clr")
  expect_identical(input, original)
  expect_identical(SummarizedExperiment::assay(tse, "counts"), as(otu_table(input), "matrix"))
  expect_identical(SummarizedExperiment::assayNames(tse), c("counts", "relative", "clr"))
  info <- S4Vectors::metadata(tse)$phyloseqV2$transformations
  expect_identical(names(info), c("relative", "clr"))
  expect_identical(info$clr$source_assay, "counts")
  expect_identical(info$clr$log_base, "e")
  expect_identical(info$relative$input_scale, "counts")
  expect_identical(TreeSummarizedExperiment::rowTree(tse), phy_tree(input))
  expect_error(transform_composition(tse, "clr"), "overwrite")
  expect_error(transform_composition(tse, "relative", name="counts"), "overwrite")
  expect_error(transform_composition(tse, "clr", assay="missing"), "Assay not found")
})

test_that("relative and Hellinger transformations agree with vegan", {
  skip_if_not_installed("TreeSummarizedExperiment", "2.14.0")
  for (taxa_rows in c(TRUE, FALSE)) {
    ps <- composition_fixture(zeros=TRUE, taxa_rows=taxa_rows)
    counts <- as(otu_table(ps), "matrix")
    if (!taxa_rows) counts <- t(counts)
    for (method in c("relative", "hellinger")) {
      tse <- transform_composition(ps, method)
      expected <- vegan::decostand(t(counts),
                                   method=if (method == "relative") "total" else "hellinger")
      expected <- t(as.matrix(expected))
      attributes(expected) <- attributes(counts)
      expect_equal(SummarizedExperiment::assay(tse, method), expected, tolerance=1e-12)
      values <- SummarizedExperiment::assay(tse, method)
      expect_equal(unname(colSums(if (method == "relative") values else values ^ 2)), rep(1, 3))
    }
  }
})

test_that("CLR obeys centering, scale invariance and pairwise log ratios", {
  skip_if_not_installed("TreeSummarizedExperiment", "2.14.0")
  input <- composition_fixture()
  tse <- transform_composition(input, "clr")
  clr <- SummarizedExperiment::assay(tse, "clr")
  counts <- SummarizedExperiment::assay(tse, "counts")
  expect_true(check_composition(tse, assay="clr", scale="clr")$valid)
  expect_equal(unname(colMeans(clr)), rep(0, 3), tolerance=1e-12)
  expect_equal(unname(clr[1, ] - clr[2, ]), unname(log(counts[1, ] / counts[2, ])))
  scaled <- tse
  SummarizedExperiment::assay(scaled, "counts") <- sweep(counts, 2, c(1e-100, 7, 1e100), "*")
  scaled <- transform_composition(scaled, "clr", name="scaled_clr")
  expect_equal(SummarizedExperiment::assay(scaled, "scaled_clr"), clr, tolerance=1e-12)
  relative <- transform_composition(tse, "relative")
  relative <- transform_composition(relative, "clr", assay="relative", name="relative_clr",
                                     input_scale="proportion")
  expect_equal(SummarizedExperiment::assay(relative, "relative_clr"), clr, tolerance=1e-12)
  # For the first sample (1,2,4,8), geometric mean is sqrt(8).
  expect_equal(unname(clr[, 1]), log(2) * c(-1.5, -0.5, 0.5, 1.5))
})

test_that("CLR matches independent compositions and vegan implementations", {
  skip_if_not_installed("TreeSummarizedExperiment", "2.14.0")
  skip_if_not_installed("compositions")
  input <- composition_fixture()
  tse <- transform_composition(input, "clr")
  counts <- SummarizedExperiment::assay(tse, "counts")
  result <- SummarizedExperiment::assay(tse, "clr")
  expected <- unclass(compositions::clr(compositions::acomp(t(counts))))
  expect_equal(as.numeric(result), as.numeric(t(expected)), tolerance=1e-12)
  expected <- vegan::decostand(t(counts), method="clr", MARGIN=1)
  expect_equal(as.numeric(result), as.numeric(t(expected)), tolerance=1e-12)
})

test_that("zero policies are explicit and multiplicative replacement preserves ratios", {
  skip_if_not_installed("TreeSummarizedExperiment", "2.14.0")
  input <- composition_fixture(zeros=TRUE)
  expect_error(transform_composition(input, "clr"), "explicit zero replacement")
  expect_error(transform_composition(input, "clr", zero_method="pseudocount"), "pseudocount")
  expect_error(transform_composition(input, "clr", pseudocount=1), "requires zero_method")
  expect_error(transform_composition(input, "clr", zero_method="multiplicative"), "delta")
  for (bad in list(0, -1, Inf, NA_real_, c(1, 2))) {
    expect_error(transform_composition(input, "clr", zero_method="pseudocount", pseudocount=bad), "pseudocount")
    expect_error(transform_composition(input, "clr", zero_method="multiplicative", delta=bad), "delta")
  }
  additive <- transform_composition(input, "clr", zero_method="pseudocount", pseudocount=1)
  clr <- SummarizedExperiment::assay(additive, "clr")
  counts <- SummarizedExperiment::assay(additive, "counts")
  expect_equal(unname(clr[1, ] - clr[2, ]), unname(log((counts[1, ] + 1) / (counts[2, ] + 1))))
  multiplicative <- transform_composition(input, "clr", zero_method="multiplicative", delta=0.01)
  clr <- SummarizedExperiment::assay(multiplicative, "clr")
  proportions <- exp(clr[, 1] - max(clr[, 1]))
  proportions <- proportions / sum(proportions)
  expect_equal(unname(proportions[1]), 0.01)
  expect_equal(unname(proportions[2:4]), 0.99 * c(2, 4, 8) / 14)
  expect_equal(unname(clr[2, ] - clr[3, ]), unname(log(counts[2, ] / counts[3, ])))
  scaled <- as_tree_experiment(input)
  SummarizedExperiment::assay(scaled, "counts") <- counts * 1e50
  scaled <- transform_composition(scaled, "clr", zero_method="multiplicative", delta=0.01)
  expect_equal(SummarizedExperiment::assay(scaled, "clr"), clr, tolerance=1e-12)
  dense_zeros <- as_tree_experiment(input)
  values <- counts
  values[1:3, 1] <- 0
  SummarizedExperiment::assay(dense_zeros, "counts") <- values
  expect_error(transform_composition(dense_zeros, "clr", zero_method="multiplicative", delta=0.4), "positive mass")
})

test_that("numerical extremes, singletons and invalid inputs have explicit behavior", {
  skip_if_not_installed("TreeSummarizedExperiment", "2.14.0")
  values <- matrix(c(1e308, 1e308), 2, dimnames=list(c("a", "b"), "s"))
  tse <- transform_composition(otu_table(values, TRUE), "relative")
  expect_equal(unname(SummarizedExperiment::assay(tse, "relative")), matrix(c(0.5, 0.5), 2))
  values[,] <- c(1e-300, 1e300)
  tse <- transform_composition(otu_table(values, TRUE), "clr")
  expect_true(all(is.finite(SummarizedExperiment::assay(tse, "clr"))))
  one <- otu_table(matrix(2, dimnames=list("taxon", "s")), TRUE)
  expect_equal(as.numeric(SummarizedExperiment::assay(transform_composition(one, "relative"), "relative")), 1)
  expect_error(transform_composition(one, "clr"), "two taxa")
  expect_error(transform_composition(one, "relative", pseudocount=1), "only to CLR")
  expect_error(aitchison_distance(tse), "two samples")
  bad <- composition_fixture()
  otu_table(bad)[, 1] <- 0
  expect_error(transform_composition(bad, "relative"), "no positive abundance")
  expect_error(transform_composition(matrix(1:4, 2), "relative"), "IDs")
})

test_that("Aitchison distances have known values and preserve sample identities", {
  values <- matrix(log(2) * c(-0.5, 0.5, -1, 1), 2,
                    dimnames=list(c("a", "b"), c("01", "s2")))
  distance <- aitchison_distance(values)
  expect_equal(as.numeric(distance), log(2) / sqrt(2))
  expect_identical(attr(distance, "Labels"), c("01", "s2"))
  expect_equal(unname(as.matrix(distance)), matrix(c(0, log(2)/sqrt(2), log(2)/sqrt(2), 0), 2))
  expect_error(aitchison_distance(values + 1), "mean zero")
})

test_that("diagnostic plots preserve IDs and display empty samples", {
  input <- composition_fixture()
  counts <- as(otu_table(input), "matrix")
  counts[, 2] <- 0
  plot <- plot_composition_qc(counts)
  expect_s3_class(plot, "ggplot")
  expect_identical(levels(plot$data$Sample), c("01", "s2", "s3"))
  expect_equal(nrow(plot$data), 9)
  expect_true(all(plot$data$Value[plot$data$Sample == "s2" & plot$data$Metric != "Zero fraction"] == 0))
  expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
  counts[1, 1] <- NA_real_
  expect_error(plot_composition_qc(counts), "finite nonnegative")
})

test_that("CLR ordination aligns metadata and matches Aitchison geometry", {
  skip_if_not_installed("TreeSummarizedExperiment", "2.14.0")
  tse <- transform_composition(composition_fixture(), "clr")
  plot <- plot_composition_ordination(tse, colour="group")
  expect_s3_class(plot, "ggplot")
  expect_identical(plot$data$Sample, colnames(tse))
  expect_identical(as.character(plot$data$Colour), c("case", "control", "case"))
  # Three samples need at most two centered PCA axes, so all distances are retained.
  coordinates <- as.matrix(plot$data[c("Axis1", "Axis2")])
  expect_equal(as.numeric(stats::dist(coordinates)), as.numeric(aitchison_distance(tse)), tolerance=1e-12)
  expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
  expect_error(plot_composition_ordination(tse, colour="absent"), "metadata column")
  values <- SummarizedExperiment::assay(tse, "clr")
  expect_error(plot_composition_ordination(values, colour="group"), "metadata are required")
  values[,] <- 0
  expect_error(plot_composition_ordination(values), "variation between samples")
})

test_that("provenance detects stale or modified assays while allowing sample subsets", {
  skip_if_not_installed("TreeSummarizedExperiment", "2.14.0")
  tse <- transform_composition(composition_fixture(), "clr")
  expect_true(check_composition(tse[, c(3, 1)], assay="clr", scale="clr")$valid)
  expect_true(check_composition(tse[4:1, ], assay="clr", scale="clr")$valid)
  expect_error(aitchison_distance(tse[1:3, ]), "taxon set changed")
  changed <- tse
  SummarizedExperiment::assay(changed, "counts")[1, 1] <- 10
  expect_error(aitchison_distance(changed), "values changed")
  changed <- tse
  SummarizedExperiment::assay(changed, "clr")[1, 1] <-
    SummarizedExperiment::assay(changed, "clr")[1, 1] + 0.5
  SummarizedExperiment::assay(changed, "clr")[2, 1] <-
    SummarizedExperiment::assay(changed, "clr")[2, 1] - 0.5
  expect_error(aitchison_distance(changed), "values changed")
  changed <- tse
  SummarizedExperiment::assay(changed, "counts") <- NULL
  expect_error(aitchison_distance(changed), "source assay is missing")
  chain <- transform_composition(tse, "relative")
  chain <- transform_composition(chain, "clr", assay="relative", name="relative_clr")
  SummarizedExperiment::assay(chain, "counts")[1, 1] <- 10
  expect_error(aitchison_distance(chain, assay="relative_clr"), "values changed")
  expect_true(check_composition(SummarizedExperiment::assay(tse, "clr"), scale="clr")$valid)
})

test_that("plain SummarizedExperiment objects retain unrelated assays and metadata", {
  skip_if_not_installed("SummarizedExperiment")
  values <- matrix(c(1, 2, 3, 5, 7, 11), 3,
                    dimnames=list(c("a", "b", "c"), c("s1", "s2")))
  input <- SummarizedExperiment::SummarizedExperiment(
    assays=list(counts=values, other=values * 2), metadata=list(project="test"))
  original <- input
  result <- transform_composition(input, "clr")
  expect_identical(class(result), class(input))
  expect_identical(input, original)
  expect_identical(SummarizedExperiment::assay(result, "other"), values * 2)
  expect_identical(S4Vectors::metadata(result)$project, "test")
  expect_true(check_composition(result, assay="clr", scale="clr")$valid)
})

test_that("legacy fingerprints remain verifiable without weakening identifier checks", {
  skip_if_not_installed("TreeSummarizedExperiment", "2.14.0")
  tse <- transform_composition(composition_fixture(), "clr")
  expect_identical(S4Vectors::metadata(tse)$phyloseqV2$transformations$clr$fingerprint_version, 2L)
  legacy_hash <- function(values) {
    stats::setNames(vapply(seq_len(ncol(values)), function(i) {
      digest::digest(values[, i], algo="xxhash64")
    }, character(1)), colnames(values))
  }
  info <- S4Vectors::metadata(tse)
  info$phyloseqV2$transformations$clr$fingerprint_version <- NULL
  info$phyloseqV2$transformations$clr$source_hashes <-
    legacy_hash(SummarizedExperiment::assay(tse, "counts"))
  info$phyloseqV2$transformations$clr$result_hashes <-
    legacy_hash(SummarizedExperiment::assay(tse, "clr"))
  S4Vectors::metadata(tse) <- info
  expect_true(check_composition(tse, assay="clr", scale="clr")$valid)
  expect_true(check_composition(tse[4:1, c(3, 1)], assay="clr", scale="clr")$valid)
  SummarizedExperiment::assay(tse, "counts")[1, 1] <- 99
  expect_error(aitchison_distance(tse), "values changed")
  relative <- transform_composition(composition_fixture(), "relative")
  info <- S4Vectors::metadata(relative)
  record <- info$phyloseqV2$transformations$relative
  record$fingerprint_version <- NULL
  record$source_hashes <- legacy_hash(SummarizedExperiment::assay(relative, "counts"))
  record$result_hashes <- legacy_hash(SummarizedExperiment::assay(relative, "relative"))
  info$phyloseqV2$transformations$relative <- record
  S4Vectors::metadata(relative) <- info
  mixed <- transform_composition(relative, "clr", assay="relative", name="mixed_clr")
  expect_true(check_composition(mixed, assay="mixed_clr", scale="clr")$valid)
  SummarizedExperiment::assay(mixed, "counts")[1, 1] <- 99
  expect_error(aitchison_distance(mixed, assay="mixed_clr"), "values changed")
  fresh <- transform_composition(composition_fixture(), "clr")
  rownames(fresh) <- paste0("renamed_", rownames(fresh))
  expect_error(aitchison_distance(fresh), "taxon set changed")
  fresh <- transform_composition(composition_fixture(), "clr")
  colnames(fresh) <- paste0("renamed_", colnames(fresh))
  expect_error(aitchison_distance(fresh), "values changed")
  info <- S4Vectors::metadata(fresh)
  info$phyloseqV2$transformations$clr$fingerprint_version <- 999L
  S4Vectors::metadata(fresh) <- info
  expect_error(aitchison_distance(fresh), "Unsupported.*version")
})

test_that("diagnostic counts handle integer storage and non-finite values", {
  for (storage in c("integer", "double")) {
    values <- matrix(c(0, -1, 2, 0, NA, 0), 3,
                      dimnames=list(c("a", "b", "c"), c("s1", "s2")))
    storage.mode(values) <- storage
    report <- check_composition(values)
    expect_equal(report$samples$Observed, c(1, 0))
    expect_equal(report$samples$Zeros, c(1, 2))
    expect_equal(report$samples$Negative, c(1, 0))
    expect_equal(report$samples$NonFinite, c(0, 1))
    expect_equal(report$samples$NonInteger, c(0, 0))
    expect_identical(report$zero_taxa, "a")
  }
})
