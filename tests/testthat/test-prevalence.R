prevalence_values <- function() {
  matrix(c(1, 2, 3, 4,
           0, 5, 0, 6,
           0, 0, 0, 1,
           0, 0, 0, 0,
           1, 1, 2, 0), nrow=5, byrow=TRUE,
         dimnames=list(paste0("t", 1:5), c("01", "s2", "s3", "s4")))
}

prevalence_fixture <- function(taxa_rows=TRUE) {
  values <- prevalence_values()
  counts <- otu_table(if (taxa_rows) values else t(values), taxa_rows)
  taxonomy <- tax_table(matrix(c("A", "B", "A", "B", "C"), ncol=1,
                               dimnames=list(rownames(values), "Genus")))
  samples <- sample_data(data.frame(group=factor(c("A", "B", "A", "B")),
                                    row.names=colnames(values)))
  tree <- ape::read.tree(text="((t1:1,t2:1):1,(t3:1,(t4:1,t5:1):1):1);")
  sequences <- Biostrings::DNAStringSet(setNames(c("AAA", "CCC", "GGG", "TTT", "ACG"),
                                                rownames(values)))
  phyloseq(counts, taxonomy, samples, tree, sequences)
}

test_that("prevalence reports explicit counts, denominators and strict detection", {
  values <- prevalence_values()
  original <- values
  expected <- data.frame(Taxon=rownames(values), Present=c(4, 2, 1, 0, 3),
                          Total=4L, Prevalence=c(1, 0.5, 0.25, 0, 0.75))
  expect_identical(taxa_prevalence(values), expected)
  expect_equal(taxa_prevalence(values, detection=1)$Present, c(3, 2, 0, 0, 1))
  expect_equal(taxa_prevalence(values, detection=2)$Present, c(2, 2, 0, 0, 0))
  proportions <- sweep(values, 2, colSums(values), "/")
  expect_equal(taxa_prevalence(proportions, detection=0.2)$Present,
                c(4, 2, 0, 0, 2))
  expect_identical(values, original)
  # Zero-mass samples remain in the denominator, rather than inflating prevalence.
  empty <- cbind(values, empty=0)
  expect_identical(taxa_prevalence(empty)$Total, rep(5L, 5))
  expect_equal(taxa_prevalence(empty)$Prevalence, expected$Present / 5)
  expect_true(all(taxa_prevalence(values * 0)$Prevalence == 0))
})

test_that("fraction and sample-count thresholds include exact boundaries", {
  values <- prevalence_values()
  keep <- setNames(c(TRUE, TRUE, FALSE, FALSE, TRUE), rownames(values))
  expect_identical(filter_prevalence(values, min_prevalence=0.5, prune=FALSE), keep)
  expect_identical(filter_prevalence(values, min_samples=2, prune=FALSE), keep)
  expect_identical(filter_prevalence(values, min_prevalence=0.5), values[keep, , drop=FALSE])
  expect_identical(filter_prevalence(values, min_samples=2), values[keep, , drop=FALSE])
  expect_identical(filter_prevalence(values, min_prevalence=0), values)
  expect_identical(filter_prevalence(values, min_samples=0), values)
  expect_identical(filter_prevalence(values, min_prevalence=1), values[1, , drop=FALSE])
  expect_identical(filter_prevalence(values, min_samples=4), values[1, , drop=FALSE])
  expect_identical(filter_prevalence(values, min_samples=2, detection=1, prune=FALSE),
                   setNames(c(TRUE, TRUE, FALSE, FALSE, FALSE), rownames(values)))
  # 0.07 * 100 is slightly above 7 on common R platforms: no ceiling is needed.
  boundary <- rbind(seven=c(rep(1, 7), rep(0, 93)),
                     six=c(rep(1, 6), rep(0, 94)))
  colnames(boundary) <- paste0("s", 1:100)
  expect_identical(filter_prevalence(boundary, min_prevalence=0.07, prune=FALSE),
                   c(seven=TRUE, six=FALSE))
})

test_that("singleton dimensions and all-removed selections remain explicit", {
  values <- prevalence_values()
  expect_identical(filter_prevalence(values[1, , drop=FALSE], min_samples=1),
                   values[1, , drop=FALSE])
  expect_identical(filter_prevalence(values[, 1, drop=FALSE], min_samples=1),
                   values[c(1, 5), 1, drop=FALSE])
  expect_identical(filter_prevalence(values[4, , drop=FALSE], min_prevalence=0),
                   values[4, , drop=FALSE])
  expect_identical(filter_prevalence(values, min_samples=1, detection=10, prune=FALSE),
                   setNames(rep(FALSE, 5), rownames(values)))
  expect_error(filter_prevalence(values, min_samples=1, detection=10), "No taxa.*prune=FALSE")
  # Removing taxa can create empty samples; the filter never removes samples.
  filtered <- filter_prevalence(values, min_samples=2, detection=4)
  expect_identical(colnames(filtered), colnames(values))
  expect_identical(rownames(filtered), "t2")
  expect_false(check_composition(filtered)$valid)
})

test_that("phyloseq filtering aligns all components in both table orientations", {
  for (taxa_rows in c(TRUE, FALSE)) {
    input <- prevalence_fixture(taxa_rows)
    original <- input
    selected <- c("t1", "t2", "t5")
    filtered <- filter_prevalence(input, min_prevalence=0.5)
    expect_s4_class(filtered, "phyloseq")
    expect_true(methods::validObject(filtered))
    expect_identical(taxa_are_rows(filtered), taxa_rows)
    expect_identical(taxa_names(filtered), selected)
    expect_identical(sample_data(filtered), sample_data(input))
    expect_identical(tax_table(filtered), tax_table(input)[selected, , drop=FALSE])
    expect_identical(phy_tree(filtered), ape::drop.tip(phy_tree(input), c("t3", "t4")))
    expect_identical(refseq(filtered), refseq(input)[selected])
    expected <- prevalence_values()[selected, , drop=FALSE]
    expect_identical(as(otu_table(filtered), "matrix"), if (taxa_rows) expected else t(expected))
    expect_identical(taxa_prevalence(input), taxa_prevalence(prevalence_values()))
    expect_identical(filter_prevalence(input, min_samples=0), input)
    table <- otu_table(input)
    filtered_table <- filter_prevalence(table, min_samples=2)
    expect_s4_class(filtered_table, "otu_table")
    expect_identical(taxa_are_rows(filtered_table), taxa_rows)
    expect_identical(as(filtered_table, "matrix"), if (taxa_rows) expected else t(expected))
    singleton_table <- filter_prevalence(table, min_prevalence=1)
    expect_identical(dim(singleton_table), if (taxa_rows) c(1L, 4L) else c(4L, 1L))
    expect_identical(taxa_names(singleton_table), "t1")
    expect_identical(input, original)
  }
  # A single retained phyloseq taxon follows the existing documented tree rule.
  expect_warning(singleton <- filter_prevalence(prevalence_fixture(), min_prevalence=1),
                  "tree replaced with NULL")
  expect_null(phy_tree(singleton, errorIfNULL=FALSE))
  expect_true(methods::validObject(singleton))
})

test_that("TSE and SE filtering preserves assays and exposes stale transformations", {
  skip_if_not_installed("TreeSummarizedExperiment", "2.14.0")
  tse <- as_tree_experiment(prevalence_fixture())
  tse <- transform_composition(tse, "relative")
  tse <- transform_composition(tse, "clr", zero_method="pseudocount", pseudocount=1)
  original <- tse
  filtered <- filter_prevalence(tse, min_samples=2)
  selected <- c("t1", "t2", "t5")
  expect_s4_class(filtered, "TreeSummarizedExperiment")
  expect_true(methods::validObject(filtered))
  expect_identical(rownames(filtered), selected)
  expect_identical(TreeSummarizedExperiment::rowLinks(filtered)$nodeLab, selected)
  expect_identical(TreeSummarizedExperiment::referenceSeq(filtered),
                   TreeSummarizedExperiment::referenceSeq(tse)[selected])
  expect_identical(TreeSummarizedExperiment::rowTree(filtered),
                   TreeSummarizedExperiment::rowTree(tse))
  expect_identical(SummarizedExperiment::colData(filtered), SummarizedExperiment::colData(tse))
  expect_identical(SummarizedExperiment::rowData(filtered),
                   SummarizedExperiment::rowData(tse)[selected, , drop=FALSE])
  expect_identical(SummarizedExperiment::assayNames(filtered), SummarizedExperiment::assayNames(tse))
  for (assay in SummarizedExperiment::assayNames(tse)) {
    expect_identical(SummarizedExperiment::assay(filtered, assay),
                     SummarizedExperiment::assay(tse, assay)[selected, , drop=FALSE])
  }
  expect_identical(S4Vectors::metadata(filtered), S4Vectors::metadata(tse))
  expect_true(check_composition(filtered, scale="counts")$valid)
  for (assay in c("relative", "clr")) {
    report <- check_composition(filtered, assay=assay, scale=if (assay == "clr") "clr" else "proportion")
    expect_false(report$valid)
    expect_match(paste(report$issues, collapse=" "), "taxon set changed")
  }
  expect_error(taxa_prevalence(filtered, assay="relative"), "taxon set changed")
  recomputed <- transform_composition(filtered, "clr", name="filtered_clr",
                                      zero_method="pseudocount", pseudocount=1)
  expect_true(check_composition(recomputed, assay="filtered_clr", scale="clr")$valid)
  expect_identical(filter_prevalence(tse, min_samples=2, assay="relative"), filtered)
  expect_identical(tse, original)
  standard <- SummarizedExperiment::SummarizedExperiment(
    assays=list(counts=prevalence_values(), alternative=prevalence_values() * 2))
  result <- filter_prevalence(standard, min_samples=2)
  expect_identical(class(result), class(standard))
  expect_identical(result, standard[selected, , drop=FALSE])
})

test_that("sample and taxon permutations retain prevalence and alignment", {
  values <- prevalence_values()
  permuted <- values[c(5, 3, 1, 4, 2), c(4, 1, 3, 2), drop=FALSE]
  expect_equal(taxa_prevalence(permuted)$Prevalence, c(0.75, 0.25, 1, 0, 0.5))
  expect_identical(filter_prevalence(permuted, min_samples=2),
                   permuted[c(1, 3, 5), , drop=FALSE])
  skip_if_not_installed("TreeSummarizedExperiment", "2.14.0")
  tse <- transform_composition(prevalence_fixture(), "relative")
  tse <- tse[c(5, 3, 1, 4, 2), c(4, 1, 3, 2)]
  expect_identical(taxa_prevalence(tse, assay="relative"), taxa_prevalence(permuted))
})

test_that("invalid thresholds, input values and identifiers are rejected", {
  values <- prevalence_values()
  expect_error(filter_prevalence(values), "exactly one")
  expect_error(filter_prevalence(values, min_prevalence=0.5, min_samples=2), "exactly one")
  bad_numeric <- list(NA_real_, NaN, Inf, -Inf, -1, numeric(), c(0, 1), "1", TRUE, 1i)
  for (bad in bad_numeric) {
    expect_error(taxa_prevalence(values, detection=bad), "detection")
    expect_error(filter_prevalence(values, min_prevalence=bad), "min_prevalence")
    expect_error(filter_prevalence(values, min_samples=bad), "min_samples")
  }
  expect_error(filter_prevalence(values, min_prevalence=1.1), "min_prevalence")
  expect_error(filter_prevalence(values, min_samples=1.5), "whole number")
  expect_error(filter_prevalence(values, min_samples=5), "number of samples")
  for (bad in list(NA, logical(), c(TRUE, FALSE), 1, "yes")) {
    expect_error(filter_prevalence(values, min_samples=1, prune=bad), "prune")
  }
  for (bad in c(NA_real_, NaN, Inf, -Inf, -1)) {
    current <- values
    current[1, 1] <- bad
    expect_error(taxa_prevalence(current), "finite, nonnegative")
  }
  expect_error(taxa_prevalence(values + 1i), "real numeric")
  expect_error(taxa_prevalence(matrix("1", 1, 1)), "real numeric")
  expect_error(taxa_prevalence(matrix(numeric(), 0, 1)), "non-empty")
  expect_error(taxa_prevalence(unname(values)), "taxon IDs")
  for (ids in list(c("t1", "t1", "t3", "t4", "t5"), c("", "t2", "t3", "t4", "t5"),
                  c(NA, "t2", "t3", "t4", "t5"))) {
    current <- values
    rownames(current) <- ids
    expect_error(taxa_prevalence(current), "taxon IDs")
  }
  current <- values
  colnames(current)[1] <- colnames(current)[2]
  expect_error(taxa_prevalence(current), "sample IDs")
})

test_that("recorded CLR and invalidated source assays cannot define prevalence", {
  skip_if_not_installed("TreeSummarizedExperiment", "2.14.0")
  values <- prevalence_values()
  values[,] <- 1
  tse <- transform_composition(otu_table(values, TRUE), "clr")
  # A valid all-zero CLR passes nonnegativity but must never imply absence.
  expect_true(all(SummarizedExperiment::assay(tse, "clr") == 0))
  expect_error(taxa_prevalence(tse, assay="clr"), "CLR coordinates")
  expect_error(filter_prevalence(tse, min_samples=1, assay="clr"), "CLR coordinates")
  expect_error(taxa_prevalence(tse, assay="absent"), "Assay not found")
  tse <- transform_composition(tse, "relative")
  SummarizedExperiment::assay(tse, "counts")[1, 1] <- 2
  expect_error(taxa_prevalence(tse, assay="relative"), "values changed")
})
