# Compositional workflow in phyloseqV2

The new functions keep original data and derived assays in a
`TreeSummarizedExperiment`. The existing R package name remains `phyloseq`.
Conversion also retains taxonomy, sample metadata, phylogenetic tree and reference
sequences. It refuses inconsistent component IDs instead of dropping observations.

## Quick start

Use the local fork in a fresh R session:

```r
pkgload::load_all(".")
data(GlobalPatterns)

# Filter the original counts first. This 10% threshold is illustrative.
head(taxa_prevalence(GlobalPatterns))
ps <- filter_prevalence(GlobalPatterns, min_prevalence=0.10)
report <- check_composition(ps, scale="counts")
stopifnot(report$valid)
head(report$samples)

tse <- as_tree_experiment(ps)
tse <- transform_composition(tse, "relative", input_scale="counts")
tse <- transform_composition(tse, "hellinger", input_scale="counts")

# Example parameter, not a universally recommended pseudocount:
tse <- transform_composition(tse, "clr", zero_method="pseudocount",
                              pseudocount=1, input_scale="counts")
SummarizedExperiment::assayNames(tse)
# counts, relative, hellinger, clr

check_composition(tse, assay="clr", scale="clr", stop_on_error=TRUE)
d <- aitchison_distance(tse)
plot_composition_qc(tse)
plot_composition_ordination(tse, colour="SampleType")
```

`transform_composition(ps, ...)` also accepts phyloseq directly and returns a
TreeSummarizedExperiment. Ordinary SummarizedExperiment subclasses are supported
without replacing their class or unrelated assays. Counts are never overwritten;
an existing output assay name causes an error. Choose a new name for comparisons.

## Prevalence filtering

`taxa_prevalence()` reports each taxon's ID, number of detected samples (`Present`),
total number of samples (`Total`) and fraction detected (`Prevalence`).
`filter_prevalence()` keeps taxa meeting an explicit minimum fraction **or**
number of samples. Supply exactly one; no retention threshold is chosen for you.

```r
# Positive abundance in at least 10% of all samples:
ps <- filter_prevalence(GlobalPatterns, min_prevalence=0.10)

# Alternatively, more than 5 reads in at least 3 samples:
ps <- filter_prevalence(GlobalPatterns, min_samples=3, detection=5)

# Inspect a selection without modifying or pruning the input:
keep <- filter_prevalence(GlobalPatterns, min_prevalence=0.10, prune=FALSE)
table(keep)

# The same API accepts TreeSummarizedExperiment and SummarizedExperiment:
counts_tse <- as_tree_experiment(GlobalPatterns)
filtered_tse <- filter_prevalence(counts_tse, min_prevalence=0.10, assay="counts")
```

These are example thresholds, not universal recommendations. Presence requires
**abundance > detection** (strict); retention requires **prevalence >= minimum**
(inclusive). Thus `detection=0` counts any positive value and `min_prevalence=1`
requires detection in every sample. Detection uses the input assay's units:
read-count and relative-abundance thresholds are not interchangeable.

All samples, including empty samples, contribute to the denominator and remain
in the result. Missing/non-finite or negative abundances and ambiguous IDs cause
an error. If no taxa pass, pruning errors; `prune=FALSE` can still return an
all-FALSE selection. Filtering can create empty samples; inspect them with
`check_composition()` before transformation. No renormalization is implicit.

phyloseq inputs retain their table orientation and matching taxonomy, metadata
and reference sequences. Trees use native `prune_taxa()` behavior: reducing a
phyloseq tree to one tip replaces it with NULL with a warning. Experiment inputs
retain all assays and native tree/row links through row subsetting.

Filter original abundances **before zero replacement** and compositional
transformation. A pseudocount can make previously absent taxa positive. Recorded
CLR assays are rejected as a prevalence source, even if all their coordinates
happen to be zero. Nonnegative external assays without provenance cannot reveal
whether zeros were replaced; the caller must establish their suitability.
When filtering an experiment that already has transformations, records remain
unchanged so the checks below detect the changed taxon set and request
recomputation. The filter does not automatically recalculate other assays.

## What each representation means

| Assay | Definition per sample | Intended interpretation |
| --- | --- | --- |
| counts | Original observed counts | Sequencing observations, not absolute microbial load |
| relative | Counts divided by their total | Proportions summing to one; closure alone does not remove compositional effects |
| hellinger | Square root of proportions | Ecological transformation with unit squared norm; not a log-ratio method |
| clr | Natural log of each positive component minus the mean log | Relative to the geometric mean of the retained components; negative coordinates are valid |

The QC plot shows sample totals, positive entries and zero fraction on the
selected nonnegative assay. Empty samples remain visible. It does not melt the
entire abundance table or convert numeric-looking IDs such as `"01"` to numbers.
The ordination plot uses centered, unscaled PCA on CLR coordinates, with equal
axis units; axis labels show explained variance. Distances between full CLR vectors are Aitchison
distances. Two plotted axes are a projection and need not retain all distances.

The classical covariance matrix of CLR coordinates is singular because of the
zero-sum constraint. These transformations do not by themselves validate a
differential abundance model, remove batch effects, correct contamination, or
recover absolute abundances. ILR/PhILR and model-based normalization require
separate APIs and validation.

## Zeros require an explicit assumption

CLR rejects zeros by default. The function never chooses a pseudocount or a
replacement magnitude automatically, and never treats negative data as counts.

**Additive pseudocount.** Adds the chosen positive value to every component on
the input scale, then computes CLR. Adding 1 to counts is different from adding
1 to proportions. Fixed additive pseudocounts break scale invariance and alter
ratios among observed components.

**Multiplicative replacement.** Closes each sample to proportions, replaces each
zero by an explicit `delta`, and multiplies positive components by
`1 - number_of_zeros * delta`. This preserves ratios among the positive
components, but the choice of delta still changes distances involving zeros.

```r
tse <- transform_composition(tse, "clr", name="clr_delta_001",
                              zero_method="multiplicative", delta=0.001)
d_delta <- aitchison_distance(tse, assay="clr_delta_001")
```

`delta` is a proportion, not a read count, and must leave positive mass for the
observed components. Neither method can determine which zeros represent real
biological absences, sampling zeros or detection limits. Compare scientifically
defensible zero treatments as a sensitivity analysis. Avoid transformations
driven solely by obtaining a visually attractive separation.

## What the checks enforce

`check_composition()` reports issues without changing data. Use
`stop_on_error=TRUE` to stop a pipeline. It checks:

- Unique, nonempty, nonmissing sample and taxon IDs.
- Finite real values; nonnegativity and positive sample mass for abundances.
- Whole-number values when `scale="counts"` is requested.
- Unit sample sums for proportions and zero sample means for CLR coordinates.
- All-zero taxa, reported without removal.
- For transformations created here: unchanged source/output values and taxon
  sets, using recorded per-sample `digest` fingerprints, including source chains.

Numerical checks cannot prove a transformation's biological appropriateness.
External CLR assays with no record receive numerical checks only. Fingerprints
detect accidental changes; they are not a security or scientific provenance
standard. The per-sample `Valid` column describes numerical requirements; consult
the report-level `valid` and `issues` for identifier or provenance failures.

Sample subsets and sample/taxon permutations remain valid. Changing the taxon
set after CLR changes the reference composition, so the old assay is rejected
even if it happens to remain centered. Changed source counts also invalidate
derived assays. Recompute from the desired source with a new output name:

```r
filtered <- tse[taxa_to_keep, ]
check_composition(filtered, assay="clr", scale="clr")  # reports stale CLR
filtered <- transform_composition(filtered, "clr", name="clr_filtered",
                                   zero_method="pseudocount", pseudocount=1)
```

Record parameters and taxon/sample IDs are available at:

```r
S4Vectors::metadata(tse)$phyloseqV2$transformations
```

## Dependencies and performance

See the [plotting integration guide](plot-integrations.md) for optional microViz,
tidyplots and pheatmap backends, annotation exports and display-only CLR subsets.

`digest` is a small direct dependency for change detection.
`TreeSummarizedExperiment`, `SummarizedExperiment` and `S4Vectors` are optional
dependencies for the new containers; they are already installed in the current
development environment. `compositions` is used as an independent test reference.
CI installs the optional packages so these integration tests execute.

The numerical kernels use vectorized R operations and compiled base-R routines;
there is no new C++/Rust/Python build requirement. Profile representative data
before adding a compiled kernel. Assays are currently materialized as dense
matrices: using a sparse, delayed or disk-backed input does not guarantee bounded
memory use. Block processing is a separate future optimization.

After profiling, fingerprint serialization and diagnostic allocations were
optimized. A paired local benchmark on R 4.4.3 (macOS arm64) for 20,000 taxa by
100 samples reduced median CLR time from 0.258 to 0.125 seconds, including
validation and fingerprint recording. See the
[change list and full measurements](changes-and-performance.md) for methodology,
allocation measurements, limits and reproduction commands. This is not a
guarantee for larger data, other hardware, or sparse/delayed assays.

New transformations use fingerprint version 2, hashing numeric columns without
serializing taxon names repeatedly. ID sets and canonical ordering are checked
separately. Version-1 records remain verifiable with their original algorithm;
they are not rewritten implicitly.

The result uses standard TreeSummarizedExperiment assays and metadata and can be
used with `mia`/`miaViz` when installed. This change does not add mia-specific
wrappers or claim that every external operation preserves the recorded provenance.

## Mathematical and integration tests

```sh
Rscript --vanilla -e 'testthat::test_local(".", filter="composition", reporter="summary", stop_on_failure=TRUE)'
Rscript --vanilla -e 'testthat::test_local(".", filter="prevalence", reporter="summary", stop_on_failure=TRUE)'
```

Tests compare CLR against `compositions` and `vegan`, and relative/Hellinger
against `vegan`. Additional properties cover known numerical answers, centering,
scale invariance, log-ratio preservation, zero policies, extreme magnitudes,
singletons, metadata alignment, immutability of source data, Aitchison distances,
PCA geometry and stale-assay detection. Passing these tests establishes numerical
and interface behavior, not a benchmark of downstream statistical power.

Prevalence tests cover strict detection and inclusive retention boundaries,
fraction/count thresholds, all-zero samples in the denominator, invalid values
and IDs, singleton dimensions, both phyloseq orientations, preservation of
metadata/tree links/reference sequences, source immutability and stale-assay
detection after filtering. Recorded all-zero CLR assays are explicitly rejected
as a source of presence/absence information.

References: [CLR documentation](https://search.r-project.org/CRAN/refmans/compositions/html/clr.html),
[vegan transformations](https://vegandevs.github.io/vegan/reference/decostand.html),
[multiplicative replacement](https://doi.org/10.1023/A:1023866030544),
[TreeSummarizedExperiment](https://bioconductor.org/packages/TreeSummarizedExperiment/),
[mia converters and assay conventions](https://microbiome.github.io/mia/reference/convert.html).
