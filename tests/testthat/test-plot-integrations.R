plot_integration_values <- function() {
  matrix(c(1, 2, 4, 8, 2, 0, 1, 3, 6, 3, 2, 1), nrow=4,
         dimnames=list(c("01", "t2", "t3", "t4"), c("01", "s2", "s3")))
}

plot_integration_fixture <- function(taxa_rows=TRUE) {
  values <- plot_integration_values()
  counts <- otu_table(if (taxa_rows) values else t(values), taxa_rows)
  samples <- sample_data(data.frame(group=factor(c("case", "control", "case")),
                                    age=c(30, 40, 50), detected=c(TRUE, FALSE, TRUE),
                                    row.names=colnames(values)))
  taxonomy <- tax_table(matrix(c("A", "A", "B", "C"), ncol=1,
                               dimnames=list(rownames(values), "Genus")))
  tree <- ape::read.tree(text="(('01':1,t2:1):1,(t3:1,t4:1):1);")
  sequences <- Biostrings::DNAStringSet(setNames(c("AAA", "CCC", "GGG", "TTT"),
                                                rownames(values)))
  phyloseq(counts, taxonomy, samples, tree, sequences)
}

test_that("plot data retain zeros, numeric-looking IDs, values and aligned annotations", {
  for (taxa_rows in c(TRUE, FALSE)) {
    input <- plot_integration_fixture(taxa_rows)
    original <- input
    data <- composition_data(input, sample_vars=c("group", "age"), tax_vars="Genus")
    expect_identical(data$Sample, rep(c("01", "s2", "s3"), each=4))
    expect_identical(data$Taxon, rep(c("01", "t2", "t3", "t4"), times=3))
    expect_identical(data$Abundance, as.vector(plot_integration_values()))
    expect_identical(data$Sample_group, rep(factor(c("case", "control", "case")), each=4))
    expect_identical(data$Sample_age, rep(c(30, 40, 50), each=4))
    expect_identical(data$Tax_Genus, rep(c("A", "A", "B", "C"), times=3))
    expect_identical(input, original)
    subset <- composition_data(input, taxa=c("t4", "01"))
    expect_identical(subset$Taxon, rep(c("t4", "01"), times=3))
    expect_identical(subset$Abundance, as.vector(plot_integration_values()[c("t4", "01"), ]))
  }
  input <- plot_integration_fixture()
  meta <- as(sample_data(input), "data.frame")
  meta$Sample <- c("a", "b", "c")
  meta$Abundance <- 101:103
  sample_data(input) <- sample_data(meta)
  data <- composition_data(input, sample_vars=c("Sample", "Abundance"))
  expect_identical(data$Sample, rep(c("01", "s2", "s3"), each=4))
  expect_identical(data$Sample_Sample, rep(c("a", "b", "c"), each=4))
  expect_identical(data$Sample_Abundance, rep(101:103, each=4))
})

test_that("plot export validates the complete source before selecting display taxa", {
  values <- plot_integration_values()
  for (bad in c(NA_real_, NaN, Inf, -1)) {
    invalid <- values
    invalid[4, 3] <- bad
    expect_error(composition_data(invalid, taxa="01"), "Non-finite|Negative")
  }
  expect_error(composition_data(unname(values)), "taxon IDs")
  for (bad in list(character(), c("01", "01"), NA_character_, "absent", 1)) {
    expect_error(composition_data(values, taxa=bad), "taxa")
  }
  for (bad in list(NA_character_, c("age", "age"), "", 1, "absent")) {
    expect_error(composition_data(plot_integration_fixture(), sample_vars=bad), "sample_vars")
    expect_error(composition_data(plot_integration_fixture(), tax_vars=bad), "tax_vars")
  }
  expect_error(composition_data(values, sample_vars="age"), "annotation")
  one <- values["01", "s2", drop=FALSE]
  expect_identical(composition_data(one), data.frame(Sample="s2", Taxon="01", Abundance=2))
})

test_that("experiment plot data preserve reordered metadata and full-reference CLR values", {
  skip_if_not_installed("TreeSummarizedExperiment", "2.14.0")
  input <- transform_composition(plot_integration_fixture(), "relative")
  input <- transform_composition(input, "clr", zero_method="pseudocount", pseudocount=1)
  original <- input
  input <- input[c("t4", "01", "t2", "t3"), c("s3", "01", "s2")]
  selected <- c("t2", "01")
  data <- composition_data(input, assay="clr", input_scale="clr", taxa=selected,
                           sample_vars="age", tax_vars="Genus")
  expected <- SummarizedExperiment::assay(input, "clr")[selected, , drop=FALSE]
  expect_identical(data$Abundance, as.vector(expected))
  expect_identical(data$Sample, rep(c("s3", "01", "s2"), each=2))
  expect_identical(data$Sample_age, rep(c(50, 30, 40), each=2))
  expect_identical(data$Tax_Genus, rep(c("A", "A"), times=3))
  expect_false(all(abs(colMeans(expected)) < 1e-8))
  expect_true(check_composition(input, assay="clr", scale="clr")$valid)
  stale <- original[1:3, ]
  expect_error(composition_data(stale, assay="clr", input_scale="clr"), "taxon set changed")
  SummarizedExperiment::assay(input, "counts")[1, 1] <- 99
  expect_error(composition_data(input, assay="clr", input_scale="clr"), "values changed")
})

test_that("tidyplots bars preserve selected proportions without reclosing them", {
  skip_if_not_installed("tidyplots", "0.4.0")
  values <- plot_integration_values()
  proportions <- sweep(values, 2, colSums(values), "/")
  selected <- c("t4", "01")
  p <- plot_composition_tidy(proportions, input_scale="proportion", taxa=selected,
                             width=90, height=60)
  expect_s3_class(p, "tidyplot")
  expect_identical(p$data$Abundance, as.vector(proportions[selected, , drop=FALSE]))
  expect_identical(levels(p$data$Sample), colnames(values))
  expect_identical(levels(p$data$Taxon), selected)
  built <- ggplot2::ggplot_build(p)$data[[1]]
  totals <- tapply(built$ymax - built$ymin, built$x, sum)
  expect_equal(as.numeric(totals), unname(colSums(proportions[selected, , drop=FALSE])))
  expect_true(all(totals < 1))
  expect_s3_class(tidyplots::adjust_size(p, width=110), "tidyplot")
  expect_error(plot_composition_tidy(values, y="another"), "mappings are fixed")
  expect_error(plot_composition_tidy(values, input_scale="clr"), "arg")
})

test_that("pheatmap clusters unchanged values and supports annotations and singleton displays", {
  skip_if_not_installed("pheatmap", "1.0.12")
  input <- plot_integration_fixture(FALSE)
  original <- input
  h <- plot_composition_heatmap(input, sample_vars=c("group", "age", "detected"),
                                tax_vars="Genus", silent=TRUE)
  expect_s3_class(h, "pheatmap")
  for (field in c("height", "merge", "labels", "order")) {
    expect_equal(h$tree_row[[field]], stats::hclust(stats::dist(plot_integration_values()), method="complete")[[field]])
    expect_equal(h$tree_col[[field]], stats::hclust(stats::dist(t(plot_integration_values())), method="complete")[[field]])
  }
  expect_true(all(c("col_annotation", "row_annotation", "annotation_legend") %in% h$gtable$layout$name))
  annotation <- h$gtable$grobs[[which(h$gtable$layout$name == "col_annotation")]]
  expect_identical(rownames(annotation$gp$fill), h$tree_col$labels[h$tree_col$order])
  expect_identical(input, original)
  expect_error(plot_composition_heatmap(input, annotation_col=NA), "managed by the wrapper")
  expect_error(plot_composition_heatmap(input, mat=plot_integration_values()), "managed by the wrapper")
  for (value in c(1, 9)) {
    single <- matrix(value, 1, 1, dimnames=list("t", "01"))
    one <- plot_composition_heatmap(single, silent=TRUE)
    expect_s3_class(one, "pheatmap")
    expect_true(is.na(one$tree_row))
    expect_true(is.na(one$tree_col))
  }
  constant <- matrix(2, 2, 2, dimnames=list(c("a", "b"), c("01", "02")))
  expect_s3_class(plot_composition_heatmap(constant, silent=TRUE, color=c("white", "blue")), "pheatmap")
  expect_s3_class(plot_composition_heatmap(plot_integration_values()[1, , drop=FALSE], silent=TRUE), "pheatmap")
  expect_s3_class(plot_composition_heatmap(plot_integration_values()[, 1, drop=FALSE], silent=TRUE), "pheatmap")
})

test_that("CLR heatmaps validate before selecting a display subset", {
  skip_if_not_installed("pheatmap", "1.0.12")
  skip_if_not_installed("TreeSummarizedExperiment", "2.14.0")
  input <- transform_composition(plot_integration_fixture(), "clr", zero_method="pseudocount", pseudocount=1)
  selected <- c("t4", "01")
  h <- plot_composition_heatmap(input, assay="clr", input_scale="clr", taxa=selected, silent=TRUE)
  expected <- SummarizedExperiment::assay(input, "clr")[selected, , drop=FALSE]
  for (field in c("height", "merge", "labels", "order")) {
    expect_equal(h$tree_col[[field]], stats::hclust(stats::dist(t(expected)), method="complete")[[field]])
  }
  stale <- input[1:3, ]
  expect_error(plot_composition_heatmap(stale, assay="clr", input_scale="clr", silent=TRUE), "taxon set changed")
  zero <- SummarizedExperiment::SummarizedExperiment(
    assays=list(clr=matrix(0, 2, 2, dimnames=list(c("a", "b"), c("01", "02")))))
  expect_s3_class(plot_composition_heatmap(zero, assay="clr", input_scale="clr", silent=TRUE), "pheatmap")
})

test_that("microViz adaptation preserves original phyloseq components and values", {
  skip_if_not_installed("microViz", "0.12.6")
  for (taxa_rows in c(TRUE, FALSE)) {
    input <- plot_integration_fixture(taxa_rows)
    original <- input
    mv <- as_microviz(input, input_scale="counts")
    expect_s4_class(mv, "psExtra")
    ps <- microViz::ps_get(mv)
    expect_identical(otu_table(ps), otu_table(input))
    expect_identical(sample_data(ps), sample_data(input))
    expect_identical(phy_tree(ps), phy_tree(input))
    expect_identical(refseq(ps), refseq(input))
    expect_identical(tax_table(ps)[, "Genus", drop=FALSE], tax_table(input))
    expect_identical(unname(as(tax_table(ps), "matrix")[, "unique"]), taxa_names(input))
    expect_identical(input, original)
  }
  # Taxonomy-free input receives an ID rank, never an invented biological rank.
  minimal <- as_microviz(otu_table(plot_integration_values(), TRUE))
  expect_identical(rank_names(microViz::ps_get(minimal)), "unique")
  expect_identical(as(otu_table(microViz::ps_get(minimal)), "matrix"), plot_integration_values())
  # All-numeric-looking IDs must survive the downstream composition plot too.
  values <- plot_integration_values()
  dimnames(values) <- list(c("01", "02", "03", "04"), c("01", "02", "03"))
  bar <- suppressMessages(microViz::comp_barplot(as_microviz(values), tax_level="unique",
                                                n_taxa=4, sample_order="default"))
  expect_setequal(as.character(bar$data$Sample), colnames(values))
  expect_setequal(as.character(bar$data$OTU), rownames(values))
  expect_equal(as.numeric(tapply(bar$data$Abundance, bar$data$Sample, sum)), rep(1, 3))
})

test_that("microViz consumes adapters for composition plots, heatmaps and CLR PCA", {
  skip_if_not_installed("microViz", "0.12.6")
  skip_if_not_installed("TreeSummarizedExperiment", "2.14.0")
  input <- transform_composition(plot_integration_fixture(), "relative")
  input <- transform_composition(input, "clr", zero_method="pseudocount", pseudocount=1)
  original <- input
  mv <- as_microviz(input, input_scale="counts")
  ps <- microViz::ps_get(mv)
  expect_identical(sample_names(ps), colnames(input))
  expect_identical(taxa_names(ps), rownames(input))
  expect_identical(as(otu_table(ps), "matrix"), SummarizedExperiment::assay(input, "counts"))
  expect_identical(as(sample_data(ps), "data.frame"), as.data.frame(SummarizedExperiment::colData(input)))
  bar <- microViz::comp_barplot(mv, tax_level="Genus", sample_order="default", n_taxa=3)
  expect_s3_class(bar, "ggplot")
  expect_equal(as.numeric(tapply(bar$data$Abundance, bar$data$Sample, sum)), rep(1, 3))
  mv_clr <- as_microviz(input, assay="clr", input_scale="clr")
  expect_identical(as(otu_table(microViz::ps_get(mv_clr)), "matrix"), SummarizedExperiment::assay(input, "clr"))
  fit <- microViz::ord_calc(mv_clr, method="PCA", scale=FALSE, verbose=FALSE)
  expected_variance <- stats::prcomp(t(SummarizedExperiment::assay(input, "clr")), scale.=FALSE)$sdev ^ 2
  expect_equal(as.numeric(vegan::eigenvals(microViz::ord_get(fit))),
                unname(expected_variance[expected_variance > 1e-12]))
  expect_s3_class(microViz::ord_plot(fit, colour="group"), "ggplot")
  heat <- microViz::comp_heatmap(mv_clr, sample_seriation="Identity", tax_seriation="Identity")
  expect_s4_class(heat, "Heatmap")
  expect_identical(heat@matrix, SummarizedExperiment::assay(input, "clr"))
  expect_identical(input, original)
  stale <- input[1:3, ]
  expect_error(as_microviz(stale, assay="clr", input_scale="clr"), "taxon set changed")
})
