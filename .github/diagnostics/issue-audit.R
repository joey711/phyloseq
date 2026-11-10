# Diagnostic reproductions, not a passing regression suite: failures are reported
# as data so all issues can be inspected in one run. Use a fresh R process.
args <- commandArgs(trailingOnly=TRUE)
root <- if (length(args)) args[1] else "."
pkgload::load_all(root, quiet=TRUE)
results <- list()
probe <- function(id, expr) {
  warnings <- character()
  value <- tryCatch(withCallingHandlers(force(expr), warning=function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error=function(e) list(error=conditionMessage(e)))
  results[[id]] <<- list(result=value, warnings=warnings)
  cat(id, ": ", jsonlite::toJSON(results[[id]], auto_unbox=TRUE), "\n", sep="")
}
counts <- otu_table(matrix(c(2, 3, 5, 7, 11, 13), nrow=3,
                          dimnames=list(c("t1", "t2", "t3"), c("s1", "s2"))), TRUE)
taxonomy <- tax_table(matrix(c("same", "same", "other"), ncol=1,
                            dimnames=list(taxa_names(counts), "Genus")))
samplemeta <- sample_data(data.frame(group=c("case", "control"),
                                    row.names=sample_names(counts)))
ps <- phyloseq(counts, taxonomy, samplemeta)

probe("1790_partial_refseq", {
  a <- phyloseq(otu_table(matrix(10, dimnames=list("asv1", "a")), TRUE),
                sample_data(data.frame(meta1="foo", meta2="bar", row.names="a")))
  b <- phyloseq(otu_table(matrix(c(10, 4), nrow=2,
                                 dimnames=list(c("asv1", "asv2"), "b")), TRUE),
                sample_data(data.frame(meta1="bar", meta2="fu", row.names="b")))
  c <- phyloseq(otu_table(a), sample_data(a),
                Biostrings::DNAStringSet(c(asv1="ACTC")))
  plain <- merge_phyloseq(a, b)
  with_refseq <- merge_phyloseq(c, b)
  list(expected_taxa=ntaxa(plain), actual_taxa=ntaxa(with_refseq),
       expected_reads=sum(sample_sums(plain)), actual_reads=sum(sample_sums(with_refseq)))
})
probe("1744_single_metadata_column", {
  a <- phyloseq(otu_table(matrix(10, dimnames=list("asv1", "a")), TRUE),
                sample_data(data.frame(group="case", row.names="a")))
  b <- phyloseq(otu_table(matrix(20, dimnames=list("asv1", "b")), TRUE),
                sample_data(data.frame(group="control", row.names="b")))
  nsamples(merge_phyloseq(a, b))
})
probe("1359_single_taxonomy_column", ntaxa(tax_glom(ps, "Genus")))
probe("additional_single_sample_merge_taxa", {
  x <- otu_table(matrix(1:3, ncol=1, dimnames=list(c("t1", "t2", "t3"), "s1")), TRUE)
  as(merge_taxa(x, c("t1", "t2")), "matrix")
})
probe("873_numeric_sample_ids", {
  ids <- c("01", "10", "12", "13")
  x <- otu_table(matrix(seq_len(16), nrow=4, byrow=TRUE,
                        dimnames=list(ids, paste0("taxon", 1:4))), FALSE)
  sd <- sample_data(data.frame(group=c("A", "A", "B", "B"), row.names=ids))
  input <- phyloseq(x, sd)
  output <- psmelt(input)
  list(expected_samples=sample_names(input), actual_samples=sort(unique(output$Sample)),
       expected_reads=sum(sample_sums(input)), actual_reads=sum(output$Abundance))
})
probe("1316_subset_taxa_scope", {
  select_genus <- function(x, target) subset_taxa(x, Genus == target)
  taxa_names(select_genus(ps, "same"))
})
probe("487_subset_samples_scope", {
  select_group <- function(x, target) subset_samples(x, group == target)
  sample_names(select_group(ps, "case"))
})
probe("608_merge_samples_factors", {
  ids <- paste0("s", 1:4)
  x <- otu_table(matrix(1:12, nrow=3,
                        dimnames=list(c("t1", "t2", "t3"), ids)), TRUE)
  sd <- sample_data(data.frame(group=factor(c("A", "A", "B", "B")),
                               habitat=factor(c("soil", "soil", "water", "water")),
                               time=c(1, 1, 2, 2), row.names=ids))
  as(sample_data(merge_samples(phyloseq(x, sd), "group")), "data.frame")
})
probe("additional_rng_without_existing_state", {
  if (exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE)) {
    rm(list=".Random.seed", envir=.GlobalEnv)
  }
  sample_sums(rarefy_even_depth(counts, sample.size=3, rngseed=42,
                               replace=FALSE, trimOTUs=FALSE, verbose=FALSE))
})
probe("1262_multiargument_transform_already_fixed", {
  x <- transform_sample_counts(counts, function(x, offset) x + offset, offset=2)
  identical(as(x, "matrix"), as(counts, "matrix") + 2)
})
probe("1794_plot_deprecation", {
  plot_bar(ps)
  list(ggplot2=as.character(packageVersion("ggplot2")))
})
