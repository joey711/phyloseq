# Run from the repository root. The output directory is optional.
args <- commandArgs(trailingOnly=TRUE)
output <- if (length(args)) args[1] else file.path(tempdir(), "phyloseqV2-plots")
dir.create(output, recursive=TRUE, showWarnings=FALSE)
pkgload::load_all(".", quiet=TRUE, export_all=FALSE)
for (package in c("microViz", "tidyplots", "pheatmap")) {
  if (!requireNamespace(package, quietly=TRUE)) stop("Install optional package: ", package)
}
data(GlobalPatterns)
ps <- filter_prevalence(GlobalPatterns, min_prevalence=0.10)
ps <- prune_samples(sample_names(ps)[1:12], ps)
features <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:100]
ps <- prune_taxa(features, ps)
tse <- transform_composition(ps, "relative", input_scale="counts")
tse <- transform_composition(tse, "clr", zero_method="pseudocount", pseudocount=1,
                              input_scale="counts")
original <- tse

tidy <- plot_composition_tidy(tse, assay="relative", input_scale="proportion",
                               taxa=features[1:8], width=150, height=80)
tidy <- tidyplots::adjust_title(tidy, "tidyplots: selected taxa (8 of 100)", fontsize=10)
tidy <- tidy + ggplot2::labs(y="Relative abundance")
tidyplots::save_plot(tidy, filename=file.path(output, "tidyplots-relative.png"),
                     view_plot=FALSE, dpi=160)

mv <- as_microviz(tse, input_scale="counts")
bar <- microViz::comp_barplot(mv, tax_level="unique", n_taxa=8, sample_order="default") +
  ggplot2::labs(title="microViz: top taxa and Other",
                subtitle="Illustrative GlobalPatterns subset | 100 retained taxa, 12 samples")
ggplot2::ggsave(file.path(output, "microviz-composition.png"), bar,
                width=9, height=5.5, units="in", dpi=160, bg="white")

selected <- features[1:16]
values <- SummarizedExperiment::assay(tse, "clr")[selected, , drop=FALSE]
limit <- max(abs(values))
heatmap <- plot_composition_heatmap(
  tse, assay="clr", input_scale="clr", taxa=selected,
  sample_vars="SampleType", tax_vars="Phylum",
  color=grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks=seq(-limit, limit, length.out=101),
  main="pheatmap: CLR display subset | pseudocount = 1",
  filename=file.path(output, "pheatmap-clr.png"), width=9, height=6.5,
  silent=TRUE)

mv_clr <- as_microviz(tse, assay="clr", input_scale="clr")
fit <- microViz::ord_calc(mv_clr, method="PCA", scale=FALSE, verbose=FALSE)
ord <- microViz::ord_plot(fit, colour="SampleType") +
  ggplot2::coord_equal() +
  ggplot2::labs(title="microViz: PCA of validated CLR values",
                subtitle="Illustrative GlobalPatterns subset | additive pseudocount = 1")
ggplot2::ggsave(file.path(output, "microviz-clr-pca.png"), ord,
                width=9, height=5.5, units="in", dpi=160, bg="white")

stopifnot(identical(tse, original), check_composition(tse, assay="clr", scale="clr")$valid)
writeLines(capture.output(sessionInfo()), file.path(output, "session-info.txt"))
cat("Example plots saved in:", normalizePath(output), "\n")
