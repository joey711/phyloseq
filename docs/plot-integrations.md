# microViz, tidyplots and pheatmap integration

phyloseqV2 now provides four entry points for the three optional plotting
packages. They use the same value, ID and transformation-provenance checks as
the compositional workflow. The package name remains `phyloseq`.

| Package | Main use | phyloseqV2 entry point | Returned object |
| --- | --- | --- | --- |
| microViz | Taxonomic composition plots, sample ordering, heatmaps and ordinations | `as_microviz()` | `psExtra`, consumed by microViz functions |
| tidyplots | Incrementally build and style scientific plots | `plot_composition_tidy()` | Editable `tidyplot` with stacked bars |
| pheatmap | Heatmaps with row/column clustering and metadata annotations | `plot_composition_heatmap()` | Native `pheatmap` result with dendrograms and gtable |
| Any data-frame plotting tool | Custom mappings, facets and other plot types | `composition_data()` | Long data frame with IDs, values and selected annotations |

The R package is **tidyplots**, plural; `tidyplot()` is its constructor. These
are integrations with the upstream packages, rather than bundled copies of
their implementations. Minimum optional versions are microViz 0.12.6,
tidyplots 0.4.0 and pheatmap 1.0.12. The local validation used microViz 0.12.6,
tidyplots 0.4.0 and pheatmap 1.0.13. CI is configured to install these optional
backends, with microViz from its GitHub repository; CI has not been run locally.

## Dependencies and setup

The three packages are in `Suggests`, so basic phyloseqV2 use does not require
loading them. Each entry point checks its backend and gives installation guidance
if it is missing or too old. In particular, microViz already imports phyloseq;
making it a required phyloseq import would introduce a circular dependency.

```r
# Only needed when absent; these packages are already installed locally.
install.packages(c("tidyplots", "pheatmap"))

# Follow the current microViz installation guide for Bioconductor prerequisites:
# https://david-barnett.github.io/microViz/
# Install into the library containing your intended phyloseq version.
install.packages("microViz", repos=c(
  microViz="https://david-barnett.r-universe.dev",
  CRAN="https://cloud.r-project.org"
))
```

The optional packages have their own R and dependency requirements. Use a fresh
R session when switching between the installed upstream phyloseq and this fork,
because they share a package name and S4 classes. No global package installation
is performed automatically by the new functions.

Prepare a small illustrative analysis from the repository root:

```r
pkgload::load_all(".")
data(GlobalPatterns)
ps <- filter_prevalence(GlobalPatterns, min_prevalence=0.10)
# Illustrative count-based feature selection before defining the CLR reference:
features <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:100]
ps <- prune_taxa(features, ps)

tse <- transform_composition(ps, "relative", input_scale="counts")
tse <- transform_composition(tse, "clr", zero_method="pseudocount",
                              pseudocount=1, input_scale="counts")
display_taxa <- rownames(tse)[1:12]
```

The 10% threshold, feature selection and additive pseudocount are demonstration
choices, not recommendations for every study. See the
[compositional guide](compositional-workflow.md) for their assumptions.

## tidyplots and custom data exports

```r
p <- plot_composition_tidy(tse, assay="relative", input_scale="proportion",
                            taxa=display_taxa, sample_vars="SampleType",
                            width=150, height=80)
p <- tidyplots::adjust_title(p, "Selected taxa: original relative abundances")
p
```

The bars use the selected values as supplied. Selecting 12 of 100 taxa does
**not** close the displayed fractions back to one. For a top-taxa plot that
collects the remaining taxa into an Other category, use microViz below.
CLR coordinates are rejected as input to compositional stacked bars.

Use the long data export when you need different mappings or plot types:

```r
d <- composition_data(tse, assay="relative", input_scale="proportion",
                       taxa=display_taxa, sample_vars="SampleType", tax_vars="Phylum")
head(d)
# Sample, Taxon, Abundance, Sample_SampleType, Tax_Phylum

p <- tidyplots::tidyplot(d, x=Sample_SampleType, y=Abundance, color=Taxon)
p <- tidyplots::add_mean_bar(p)
```

`composition_data()` retains zeros, numeric-looking IDs such as `"01"`, and
the selected taxon/sample order. Prefixes keep metadata names from overwriting
the core columns. Annotations are joined by IDs and retain factors. Requested
annotations must be scalar numeric, character, logical or factor columns; list
and matrix columns require an explicit upstream conversion or exclusion.
The export does not aggregate taxonomy, compute inferential statistics or change
the original data. Group summaries in a plot remain descriptive summaries.

## pheatmap with sample and taxonomy annotations

```r
values <- SummarizedExperiment::assay(tse, "clr")[display_taxa, , drop=FALSE]
limit <- max(abs(values))
h <- plot_composition_heatmap(
  tse, assay="clr", input_scale="clr", taxa=display_taxa,
  sample_vars="SampleType", tax_vars="Phylum",
  color=grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks=seq(-limit, limit, length.out=101),
  silent=TRUE
)
grid::grid.newpage()
grid::grid.draw(h$gtable)

# Alternatively pass filename="heatmap.pdf" to save with pheatmap.
```

Row and column clustering default to pheatmap's Euclidean distance and complete
linkage, with no display standardization (`scale="none"`). Clustering is disabled
automatically for a singleton dimension. Constant assays receive nondegenerate
color breaks unless explicit breaks are supplied. Additional pheatmap options
such as colors, clustering methods, labels and filenames are forwarded.

Validation occurs on the **entire source assay before display selection**.
Displaying selected CLR coordinates does not change their original reference
composition or recenter them. Distances between samples on those displayed
coordinates describe that projection; they are not the full Aitchison distances.
Use `aitchison_distance(tse)` for distances using all validated CLR coordinates.
Explicit row/column standardization changes clustering geometry; it is not a
compositional normalization. Actually filtering the TSE first invalidates the old
assay and is rejected until you recompute it from the filtered source.

## microViz composition plots, heatmaps and ordinations

```r
# The default assay is counts, validated here as whole-number counts.
mv <- as_microviz(tse, input_scale="counts")
bar <- microViz::comp_barplot(mv, tax_level="unique", n_taxa=8,
                              sample_order="default")
bar

# For biological ranks, first inspect/fix ambiguous or missing taxonomy using
# microViz::tax_fix(), then choose tax_level="Phylum", "Genus", etc.

# Hand over a CLR assay already computed with explicit zero handling:
mv_clr <- as_microviz(tse, assay="clr", input_scale="clr")
ord <- microViz::ord_calc(mv_clr, method="PCA", scale=FALSE, verbose=FALSE)
microViz::ord_plot(ord, colour="SampleType") + ggplot2::coord_equal()
microViz::comp_heatmap(mv_clr, taxa=display_taxa,
                        colors=microViz::heat_palette(sym=TRUE),
                        sample_seriation="Identity", tax_seriation="Identity")
```

The adapter uses `tax_transform(trans="identity", rank="unique",
keep_counts=FALSE)`: no extra mathematical transformation, aggregation or zero
replacement is applied. The microViz identity label means *no transformation by
the adapter*; it does not relabel the selected CLR assay as raw counts. Do not
use the CLR adapter for a count-based composition barplot. The default microViz
composition barplot closes nonnegative abundances itself, so the barplot example
uses the **count** adapter. Choose ordination methods explicitly for your data.

phyloseq inputs retain their original components, with a `unique` ID rank added
by microViz. A missing taxonomy table receives an ID rank rather than invented
biological labels. For experiment inputs, the plotting adapter transfers the
selected assay, scalar sample metadata and character/factor row annotations as
taxonomy. Numeric row annotations are not interpreted as taxonomy. Other assays,
experiment metadata, trees and reference sequences remain in the unchanged
source experiment and are not transferred. This adapter is **not a lossless
round-trip conversion**.

phyloseqV2 validates source provenance on entry and verifies that adaptation did
not change values or IDs. Subsequent microViz operations use microViz's own
metadata and checks. They do not inherit phyloseqV2 fingerprints, including its
checks against stale assays. Keep the TSE as the authoritative analysis object;
recreate a plotting adapter after changes rather than copying its values back.
Some native microViz transformations choose a pseudocount automatically: use
`transform_composition()` first when you require explicit zero-policy recording.

## Verification and runnable example

```sh
Rscript --vanilla -e 'testthat::test_local(".", filter="plot-integrations", reporter="summary", stop_on_failure=TRUE)'
Rscript --vanilla .github/diagnostics/plot-integrations-demo.R /tmp/phyloseqV2-plots
```

Tests cover both phyloseq orientations, experiment metadata alignment,
numeric-looking IDs, zero retention, invalid inputs, provenance checks before
display selection, unchanged CLR values, non-reclosed stacked bars, constant
heatmaps, singleton dimensions and upstream microViz barplot/heatmap/PCA calls.
The PCA eigenvalues are compared with an independent base-R PCA.

The adapters currently materialize dense matrices. The long export has one row
per displayed taxon/sample pair; select display taxa when this would be large.
Source validation still checks the full assay. No plotting speedup or bounded
memory guarantee is claimed for sparse/delayed inputs.

Upstream documentation: [microViz overview](https://david-barnett.github.io/microViz/),
[microViz transformations](https://david-barnett.github.io/microViz/reference/tax_transform.html),
[microViz heatmaps](https://david-barnett.github.io/microViz/reference/comp_heatmap.html),
[tidyplots](https://jbengler.github.io/tidyplots/),
[pheatmap manual](https://www.stat.ethz.ch/CRAN/web/packages/pheatmap/refman/pheatmap.html).
