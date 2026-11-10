# Changes and performance — updated 2026-11-08

This report records implementation and validation completed before publication.
Commit and push operations are managed separately by the repository owner.

## Changes implemented so far

| Area | Change | Reason |
| --- | --- | --- |
| Source structure | Extracted sample transformations, thresholded ranking and transpose into `R/transform-methods.R`; updated `DESCRIPTION` load order | Separate transformation code from filtering, pruning and rarefaction |
| Generic transformations | Removed the synthetic callback invocation; pass additional arguments on real samples and validate each result | Support multi-argument/data-dependent callbacks and avoid extra side effects |
| Data identity | Preserve dimensions, orientation and taxon/sample names even when callbacks return unnamed vectors | Prevent accidental ID changes, including singleton cases |
| Threshold ranks | Forward `keep0s` correctly from `threshrankfun()` | Honor the documented zero-handling option |
| Test structure | Deterministic fixtures, numerical checks and metadata-preservation tests; simplified package test runner | Make regressions reproducible and failures easier to interpret |
| Prevalence filtering | Added `taxa_prevalence()` diagnostics and `filter_prevalence()` for phyloseq, otu_table, SummarizedExperiment/TSE and matrices | Select taxa with explicit detection and fraction/count thresholds while retaining samples and alignment; leave stale assays detectable |
| Optional integration | Keep core filtering tests independent of `genefilter`; skip only its integration test when absent | Optional dependencies should not stop unrelated tests |
| Modern container | Added `as_tree_experiment()` | Retain values, taxonomy, sample metadata, trees and sequences in TreeSummarizedExperiment; reject inconsistent IDs |
| Compositional transforms | Added `transform_composition()` for proportions, Hellinger and CLR, with explicit CLR zero policies | Store derived assays beside original data; avoid guessing pseudocounts or zero-replacement magnitudes |
| Validation | Added `check_composition()` | Check scales, real/finite values, IDs, zero-mass samples, CLR centering and closed proportions |
| Provenance | Record transformation parameters, ID sets and per-sample source/result fingerprints; verify dependency chains | Detect transformed assays made stale by filtering or data edits |
| Distance | Added `aitchison_distance()` | Compute Euclidean distances from validated CLR coordinates without implicit preprocessing |
| Plots | Added sample QC and CLR PCA plots with aligned metadata and equal axis units | Inspect sample totals/zeros and visualize CLR geometry while preserving IDs |
| Plotting integrations | Added `composition_data()`, `plot_composition_tidy()`, `plot_composition_heatmap()` and `as_microviz()` | Use tidyplots/pheatmap/microViz with aligned annotations, explicit assay selection and provenance checks before display subsetting |
| Dependencies | Added runtime `digest`; optional TreeSummarizedExperiment/SummarizedExperiment/S4Vectors and microViz/tidyplots/pheatmap; `compositions` as a test reference | Support provenance, modern containers, plotting backends and independent numerical comparisons |
| CI and packaging | Added GitHub Actions tests, package build exclusions and eighteen public exports/help pages | Make the new APIs installable and test optional integrations in CI |
| Documentation | Added development guide, compositional workflow, diagnostic examples and issue audit | Explain assumptions, reproducible use and remaining work |

The UniFrac non-binary-tree defect and partial-refseq merge losses are addressed
in the 2026-11-08 changes below. Numeric-ID losses in legacy `psmelt()`,
categorical metadata coercion, single-column merge failures and filtering scope
issues have **not** been fixed by these changes. See
[`../.github/ISSUE_AUDIT.md`](../.github/ISSUE_AUDIT.md).

## Targeted 16S reference workflow — 2026-11-08

- Added `asv_sequences()` and `select_asvs()` for observed feature IDs or exact
  sequences, preserving numeric-looking IDs, duplicate DNA with distinct IDs,
  native containers, samples and assay provenance.
- Added `silva_reference()` with explicit releases 144/138.2, optional verified
  downloads, persistent user-selected caches and paired local-file support.
  Added runtime XVector for chunked FASTA input; no database is bundled.
- Added `filter_silva_reference()` with bacterial/archaeal domain selection,
  separate length minima, sequence/alignment and ambiguity/homopolymer/vector
  criteria, missing-value policy and per-reference exclusions. Release-144
  Pintail requests fail explicitly. Sequence/report streaming avoids loading a
  full alignment; quality metadata still occupy memory.
- Added `match_asvs_to_reference()` for full-query ungapped Biostrings searches,
  exact by default, both strands and explicit one-to-many/no-hit/truncated reports.
- Added `infer_asv_tree()` using optional DECIPHER/phangorn, joint alignment and
  JC/GTR ML inference, optional gamma/bootstrap, overlap/sequence/branch checks,
  sequence fingerprints and restored caller RNG state. Outputs separate observed
  ASVs from references; original counts/sequences are unchanged.
- Fixed UniFrac descendant accumulation for polytomies/unary nodes without
  inventing edges. Added an independent branch-incidence regression reference.
- Added `merge_phyloseq(refseq_policy=...)`: incomplete sequences now stop by
  default instead of silently deleting count-table taxa. Explicit `"drop"` warns
  and discards that component; other legacy component reconciliation is unchanged.

This is a targeted joint ML backend, not fixed-backbone placement or full-gene
reconstruction. EPA-ng/jplace, indexed gapped candidate search and additional
chimera screening remain future integrations. No performance improvement or
biological accuracy is claimed for the new phylogeny workflow. A complete
release-144 download now passes publisher checksums and whole-database QC/search,
followed by targeted joint ML inference. The live test uses reference-derived
fragments, not independently observed biological ASVs; release-138.2 full-data
checks and biological benchmarking remain additional work. See the
[usage guide](asv-silva-workflow.md) and [design roadmap](asv-silva-phylogeny-design.md).

The final isolated installed package passes **1,310 expectations, zero failures,
zero test warnings and zero skips**. `R CMD check --no-manual --ignore-vignettes
--no-examples` finishes with **Status: OK** and no package-check errors, warnings
or notes. The three new installed help examples containing executable code pass
separately. The synthetic end-to-end demo also runs successfully and preserves
its input object. Full legacy examples, vignettes, the PDF manual, remote CI and
biological accuracy validation were not run. Live release-144 integration was
then run separately; see the [report](validation/silva144-live.md). Local backends: Biostrings
2.74.1, DECIPHER 3.2.0 and phangorn 2.12.1 on R 4.4.3/macOS arm64.

## Performance optimization in this iteration

Profiling the existing compositional workflow showed that serialization inside
fingerprint creation dominated the sampled execution time (about 78% in that
specific profile). Taxon IDs were being serialized repeatedly for each sample.
The optimizations address that measured cost without disabling validation:

1. Fingerprint version 2 hashes unnamed numeric columns. Taxon ID sets and
   canonical row order are checked separately, and hashes remain keyed by sample
   ID. Numeric storage types remain part of the serialized object.
2. Records with no fingerprint-version field are verified using the original
   version-1 algorithm. Unknown versions are rejected. Existing saved assays
   remain verifiable; their original fingerprints are not silently rewritten.
3. Provenance verification avoids copying matrices into an identical row order.
   If taxa have been permuted, it still restores the recorded order before hashing.
4. Diagnostics reuse zero counts, derive positive-entry counts from disjoint
   value categories, avoid copying already-finite inputs, and avoid rounding
   integer-storage matrices to determine whether their entries are whole numbers.

No new dependency or compiled extension was needed for these optimizations.
The public APIs, scientific formulas and zero policies remain the same. New
transformation metadata contains `fingerprint_version=2L`.

## Measured results

The comparison used R 4.4.3 on macOS arm64, synthetic dense double matrices with
positive values cycling through 1 to 17, two matrix shapes, and five timed calls
per version and operation. Before/after calls alternated order, after warm-up,
with garbage collection outside the timed region. Each implementation created
its own provenance records: version 1 before, version 2 after.

Assays, IDs, diagnostics and scientific transformation parameters were compared
before timing. Fingerprint bytes intentionally differ between their versioned
representations. Both old and new representations have regression tests for
stale data and ID changes.

| Operation, 20,000 taxa × 100 samples | Before, median s | After, median s | Speedup | Allocation before → after, MiB |
| --- | ---: | ---: | ---: | ---: |
| Validate counts | 0.043 | 0.026 | 1.65× | 130.7 → 61.8 |
| Relative abundance | 0.258 | 0.130 | 1.98× | 420.3 → 300.0 |
| Hellinger | 0.259 | 0.135 | 1.92× | 435.6 → 315.3 |
| CLR | 0.258 | 0.125 | 2.06× | 382.1 → 261.9 |
| Validate CLR, including provenance | 0.216 | 0.086 | 2.51× | 323.9 → 171.6 |
| CLR from relative assay, including source-chain validation | 0.424 | 0.179 | 2.37× | 575.3 → 371.7 |

The 4,000 taxa × 500 samples shape showed similar results: CLR 0.260 → 0.127 s,
provenance-aware CLR validation 0.231 → 0.090 s, and CLR from relative abundance
0.450 → 0.188 s.

Allocations are **total bytes allocated during one profiled call**, not peak
resident RAM or retained object size. Timings are local diagnostic measurements,
not universal speed guarantees. These benchmarks do not establish performance
for zero replacement, dense conversions of sparse/delayed assays, or all other
functions in the package. Input immutability was also checked.

Raw results, package versions and source checksums:
[`benchmarks/composition-speed.json`](benchmarks/composition-speed.json).

## Validation of the performance changes

The isolated installed package passed 825 test expectations with zero failures,
two pre-existing warnings from tree-reading tests, and one optional integration
test skipped because `genefilter` is unavailable. `R CMD check --no-manual
--ignore-vignettes --no-examples` completed with no errors or package-check
warnings and two notes about that missing suggested package and its Rd links.
The complete legacy vignette/example suite was not run in this check.

Diagnostic reports were also exactly identical to the pre-optimization source
on 28 boundary cases spanning integer/double storage, zeros, negatives,
non-finite values and expected scales. Tests cover both fingerprint versions,
mixed-version dependency chains, ID changes, sample subsets and taxon permutations.

## Prevalence filtering validation

The prevalence addition contributes 148 passing expectations. The isolated
installed package now passes 973 expectations with zero failures, the same two
pre-existing warnings and one optional `genefilter` skip. The package check
completed with zero errors, zero package-check warnings and the same two notes.
Manual/vignette checks and the full legacy example suite remained excluded;
both new installed help examples were executed separately and passed.

Integration coverage includes phyloseq in both table orientations, otu_table,
matrix, SummarizedExperiment and TreeSummarizedExperiment; metadata, tree links,
reference sequences, exact threshold boundaries, empty samples, source
immutability and stale transformations are checked. A GlobalPatterns smoke test
with the illustrative 10% threshold retained 12,711 of 19,216 taxa and all 26
samples, with valid count diagnostics afterward.

Prevalence uses vectorized comparisons and row sums without adding dependencies.
The speed measurements above concern the compositional workflow, not this new
filter; no before/after speedup is claimed for prevalence filtering.

## Plotting integration validation

The plotting integrations add 123 passing expectations. The isolated installed
package passes 1,096 expectations with zero failures, the same two existing
warnings and one optional `genefilter` skip. `R CMD check --no-manual
--ignore-vignettes --no-examples` completed without errors or package-check
warnings, with the same two notes about missing `genefilter` and its Rd links.
The four new installed help examples were executed separately and passed.

Tests execute microViz barplots, ComplexHeatmap heatmaps and unscaled CLR PCA,
compare PCA eigenvalues with base R, and compare pheatmap clustering with
distances from unchanged assay values. They cover numeric-looking IDs,
annotations, singleton/constant heatmaps, display-only CLR subsets, proportion
bars without implicit reclosure, and stale assays. Local backend versions were
microViz 0.12.6, tidyplots 0.4.0 and pheatmap 1.0.13. CI now includes the optional
plotting backends but has not been run remotely.

Four example figures and their session information were generated from the
diagnostic script in `docs/examples/plot-integrations/`. See the
[integration guide](plot-integrations.md) for examples, adapter scope and source
links. Source TSEs remain authoritative: the microViz adapter transfers the
selected assay and annotations, rather than all assays, trees and provenance.
No new plotting performance benchmark is claimed.

## Warning and dependency cleanup — 2026-11-04

The two previous test warnings came from deliberately passing a non-tree file
to `read_tree()`. ape warned that the Newick terminator was absent before
returning NULL. `read_tree()` now handles only that diagnostic and the analogous
empty-input diagnostic during its Newick format probe. File-access warnings and
other parser diagnostics remain visible. `errorIfNULL=TRUE` still raises an
explicit error when no tree can be read. This is production handling of expected
parser failures, rather than blanket warning suppression in the tests.

Regression coverage verifies compressed Newick imports, Nexus imports, invalid
and empty files, strict errors, direct phylo inputs and missing-file warnings.
The documented default NULL behavior is retained.

Installed genefilter 1.88.0 and its missing dependency annotate 1.84.0 into the
existing R 4.4 site library using BiocManager with `update=FALSE`. Existing
packages were not upgraded. The dependency remains optional in DESCRIPTION;
the integration test now executes on this machine instead of skipping.
Other machines without genefilter still skip that optional integration, while
the configured CI explicitly installs it.

The installed package passed **1,109 expectations, zero failures, zero warnings
and zero skips**. `R CMD check --no-manual --ignore-vignettes --no-examples`
completed with **Status: OK**, zero package-check errors/warnings/notes, without
relaxing the suggested-dependency requirement. The full legacy example suite,
manual and vignettes remain outside this check's scope.

During installation from source, genefilter's upstream `half_range_mode.cpp`
emitted a non-fatal compiler format warning for a diagnostic Rprintf call. Its
unchanged upstream package installed and loaded successfully; the filterfun /
kOverA integration test passed. This compiler diagnostic is separate from the
two R test warnings fixed here. The sandbox also printed repository-index access
messages during dependency checking; the package check itself completed OK.

## Reproduction

From the repository root, time the current implementation:

```sh
Rscript --vanilla .github/diagnostics/benchmark-composition.R
```

To compare against a saved source snapshot:

```sh
Rscript --vanilla .github/diagnostics/benchmark-composition.R /path/to/baseline-composition.R /path/to/results.json
```

The baseline used here was saved before editing at
`/private/tmp/phyloseq-speed-before/composition-methods.R`; its MD5 is recorded in
the results JSON. That temporary snapshot is not part of the package. The script
loads baseline and current functions into separate environments in the same R
process and does not install either version over another package.

## Remaining performance work

The new assay interface still materializes data as dense matrices. Block
processing of delayed/on-disk assays and the cost of legacy `tax_glom()` and
`psmelt()` are separate optimization tasks. They require their own benchmarks and
tests for output equivalence; the speedups above should not be generalized to
those functions.
