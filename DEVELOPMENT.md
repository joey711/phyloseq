# Developing phyloseqV2

This fork currently keeps the R package name `phyloseq` and its public S4 classes.
Use a development session to work on the fork without replacing an installed
Bioconductor version:

```r
pkgload::load_all(".")
```

Start a fresh R session before switching between the installed package and the
fork, since both define the same S4 classes and methods.

## Source layout

| Area | Source | Tests |
| --- | --- | --- |
| S4 definitions and validity | `R/allClasses.R`, `R/validity-methods.R` | `test-phyloseq.R` |
| Construction and access | `R/*-class.R`, `R/*Accessors.R`, `R/assignment-methods.R` | `test-phyloseq.R`, `test-subset.R` |
| Sample transformations and transpose | `R/transform-methods.R` | `test-transform.R` |
| Compositional validation, assays and conversion | `R/composition-methods.R` | `test-composition.R` |
| Prevalence diagnostics and filtering | `R/filter-prevalence.R` | `test-prevalence.R` |
| Compositional QC and CLR ordination | `R/plot-composition.R` | `test-composition.R` |
| microViz, tidyplots and pheatmap adapters | `R/plot-integrations.R` | `test-plot-integrations.R` |
| Observed ASV selection and sequence identity | `R/asv-methods.R` | `test-asv-silva.R` |
| Versioned SILVA reference retrieval and streaming QC | `R/silva-reference.R` | `test-asv-silva.R` |
| Candidate mapping and targeted joint ML phylogeny | `R/asv-phylogeny.R` | `test-asv-silva.R` |
| Filtering, pruning and rarefaction | `R/transform_filter-methods.R` | `test-subset.R`, `test-rarefy.R` |
| Import and merging | `R/IO-methods.R`, `R/merge-methods.R` | `test-IO.R`, `test-merge.R` |
| Distances and graphics | `R/distance-methods.R`, `R/plot-methods.R` | `test-distance.R`, `test-plot.R` |
| UniFrac topology and partial-refseq merge regressions | `R/distance-methods.R`, `R/merge-methods.R` | `test-phylogeny-integrity.R` |

R source files have an explicit load order in `DESCRIPTION` under `Collate`.
Update that field whenever a module is added or split. Keep existing exports,
S4 signatures and aliases unless an API change is deliberate and documented.

## Running tests

The package imports, `pkgload` and `testthat` must be installed. Some imports
come from Bioconductor. From the repository root:

```sh
Rscript --vanilla -e 'testthat::test_local(".", reporter="summary", stop_on_failure=TRUE)'
Rscript --vanilla -e 'testthat::test_local(".", filter="transform|subset|prevalence", reporter="summary", stop_on_failure=TRUE)'
```

These commands test the working tree. `devtools::test()` is also supported.
The package-level runner in `tests/testthat-phyloseq.R` runs the installed tests
during `R CMD check`.

`genefilter` is optional: core filtering tests always run, while the integration
test skips with an explicit reason when the package is absent. The GitHub Actions
workflow installs `genefilter` and runs the entire test suite on Linux. The legacy
Travis configuration remains available for upstream reference.

To run the genefilter integration locally, install it for the active R version:

```r
BiocManager::install("genefilter", ask=FALSE, update=FALSE)
```

It is installed in the current R 4.4 library (genefilter 1.88.0). The skip remains
appropriate on machines without this optional dependency; it is not used to
bypass a failing test. Invalid-format tree-import tests now verify quiet NULL
returns or explicit errors, while missing-file diagnostics remain visible.

The optional visualization backends are also installed by CI. microViz imports
phyloseq, so keep it in `Suggests` to avoid a required dependency cycle. See the
[plotting integration guide](docs/plot-integrations.md) for adapter scope,
provenance boundaries and source-assay selection.

CI installs optional DECIPHER/phangorn for the targeted phylogeny tests. Normal
tests use offline synthetic references, never a database download. See the
[16S workflow guide](docs/asv-silva-workflow.md) for supported releases, the joint
tree/placement distinction, reference QC and explicit download behavior. The
opt-in script now passes against the complete official SILVA 144 export; see
the [live report](docs/validation/silva144-live.md). These integration checks
remain separate from routine CI and do not establish biological accuracy.

Prefer small deterministic fixtures to randomly selected subsets or large
datasets. For transformations, verify both table orientations, one-taxon and
one-sample inputs, taxon/sample names, metadata, and invalid callback results.
Test observable behavior rather than private implementation details.

## Documentation and extension sequence

Keep roxygen comments with their functions and commit the corresponding `man/`
files. Normal documentation generation uses `devtools::document()`; review the
diff because a different roxygen version can regenerate unrelated pages.

The first refactoring isolates sample transformations. Subsequent changes should
be separate, reviewable increments:

1. Separate rarefaction and filtering/pruning into dedicated modules, keeping
   existing APIs and adding regression coverage before moving implementations.
2. Introduce explicit transformation and normalization APIs with documented input
   assumptions, zero handling, and output scale. Distinguish methods operating
   per sample from methods that estimate parameters across samples.
3. Add contaminant detection as an optional integration, with aligned sample
   metadata, explicit negative controls or concentration inputs, and a diagnostic
   result that can be reviewed before pruning taxa.
4. Add a workflow vignette connecting quality control, filtering, normalization,
   and downstream analysis on a small reproducible dataset.

The first compositional workflow is now implemented; see
[the usage guide](docs/compositional-workflow.md). It preserves original assays,
adds relative/Hellinger/CLR transformations, and validates numerical properties
and assay provenance. Remaining work includes model-based normalization,
additional log-ratio methods, scalable assay processing and contaminant detection.
