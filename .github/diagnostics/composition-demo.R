# Run from the repository root, optionally passing an output directory.
args <- commandArgs(trailingOnly=TRUE)
output <- if (length(args)) args[1] else file.path(tempdir(), "phyloseqV2-example")
dir.create(output, recursive=TRUE, showWarnings=FALSE)
pkgload::load_all(".", quiet=TRUE, export_all=FALSE)
data("GlobalPatterns")
ps <- prune_samples(sample_names(GlobalPatterns)[1:12], GlobalPatterns)
keep <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:500]
ps <- prune_taxa(keep, ps)
stopifnot(check_composition(ps, scale="counts")$valid)
tse <- transform_composition(ps, "relative", input_scale="counts")
tse <- transform_composition(tse, "clr", zero_method="pseudocount", pseudocount=1,
                              input_scale="counts")
stopifnot(check_composition(tse, assay="clr", scale="clr")$valid)

qc <- plot_composition_qc(tse) +
  ggplot2::labs(title="Sample quality controls",
                subtitle="GlobalPatterns example | 500 retained taxa, 12 samples")
ordination <- plot_composition_ordination(tse, colour="SampleType") +
  ggplot2::scale_colour_brewer(palette="Dark2") +
  ggplot2::labs(title="CLR ordination in Aitchison geometry",
                subtitle="Illustrative GlobalPatterns subset | additive pseudocount = 1")
ggplot2::ggsave(file.path(output, "composition-qc.png"), qc,
                width=9, height=6.5, units="in", dpi=160, bg="white")
ggplot2::ggsave(file.path(output, "composition-ordination.png"), ordination,
                width=9, height=6, units="in", dpi=160, bg="white")
cat("Example figures saved in:", normalizePath(output), "\n")
