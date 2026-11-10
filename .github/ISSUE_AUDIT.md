# Upstream issue audit — 2026-11-03

Scope: `joey711/phyloseq`, the parent of `l-gallucci/phyloseqV2`. The fork has
issues disabled. Indexed all **741 open issues** returned by the GitHub API,
excluding pull requests; reviewed titles, selected bodies and discussion threads,
then reproduced the most consequential candidates against the local working tree.
This is a triage of the backlog, not a claim that every issue has been reproduced.

Local baseline: upstream commit `8a6c235`, package version 1.55.2, plus the local
transformation refactoring. R 4.4.3; ape 5.8-1; Biostrings 2.74.1;
ggplot2 4.0.3; Matrix 1.7-6; reshape2 1.4.5.

Priority reflects impact on results and data integrity, followed by crashes and
workflow failures. An open issue alone is not evidence that its reported bug is
still present.

## 1. UniFrac assumes binary trees and can return incorrect distances

**Priority: high. Confirmed locally; fixed in the working tree 2026-11-08.**

Sources: [#956](https://github.com/joey711/phyloseq/issues/956),
[#1574](https://github.com/joey711/phyloseq/issues/1574); related closed discussion
[#936](https://github.com/joey711/phyloseq/issues/936).

`fastUniFrac()` constructs a matrix with exactly two children per internal node.
It does not validate that the tree is binary. A rooted tree with a polytomy can
produce a matrix-recycling warning and incorrect numerical results.

Local reproduction: tree `((a:1,b:1,c:1):1,d:2);`, with sample `u` containing
only `a` and sample `v` containing only `c`.

| Distance | Expected | Before fix | After resolving with zero-length edges |
| --- | ---: | ---: | ---: |
| Unweighted | 2/3 | 1 | 2/3 |
| Weighted, normalized | 0.5 | 0.75 | 0.5 |
| Weighted, unnormalized | 2 | 3 | 2 |

The expected values follow directly from the branch lengths: the shared branch
has length 1, each unique terminal branch has length 1, and each occupied tip is
at distance 2 from the root. `ape::multi2di(random=FALSE)` adds zero-length edges;
it does not change these branch-length totals.

This confirms the non-binary-tree failure mode. It does not establish that every
QIIME/phyloseq discrepancy in #956 has the same cause; normalization defaults also
need to be aligned when comparing implementations.

Fix: accumulate child contributions along postorder edges, supporting arbitrary
node degree and unary nodes without resolving or changing biological topology.
Regression tests compare against an independent root-path branch-incidence
calculation for binary, polytomous and unary trees, both OTU orientations and
all weighted/normalization options. The values above now match expectations.
The wrapper's existing rooting behavior is outside this fix.

Location: `R/distance-methods.R`, `fastUniFrac()`.

## 2. Merging partial reference sequences silently removes taxa

**Priority: high. Confirmed locally; protected in the working tree 2026-11-08.**

Source: [#1790](https://github.com/joey711/phyloseq/issues/1790).

When one object has `refseq` and the other does not, the merged reference-sequence
component covers only part of the union of taxa. Reconstruction with `phyloseq()`
intersects all component taxon names, dropping taxa that have no sequence.

Local reproduction: merging without `refseq` retains **2 taxa and 24 reads**;
adding a sequence for only `asv1` leaves **1 taxon and 20 reads**, with no warning.

Fix: `merge_phyloseq(..., refseq_policy="error")` now stops before reconstruction
when merged sequences do not cover the merged count-table taxa. The explicit
`refseq_policy="drop"` alternative warns and removes the incomplete refseq
component, preserving the count union subject to other components' existing
reconciliation. Tests retain 2 taxa/24 reads and exercise both table orientations,
default errors and complete-reference merges. Unknown bases are never invented.
This policy change covers refseq, not incomplete trees/taxonomy generally.

Locations: `R/merge-methods.R`, `merge_phyloseq()`; `R/phyloseq-class.R`, constructor.

## 3. Numeric-looking IDs can disappear from melted data and plots

**Priority: high. Confirmed locally.**

Source: [#873](https://github.com/joey711/phyloseq/issues/873).

`reshape2::melt()` converts numeric-looking dimnames before `psmelt()` converts
them back to character. Leading zeros have already been lost at that point.
The subsequent metadata join drops the unmatched sample.

Local reproduction: samples `01`, `10`, `12`, `13` become only `10`, `12`, `13`;
the abundance total drops from **136 to 126** with no warning. This affects the
data returned by `psmelt()` and plotting functions that consume it.

Next change: preserve identifiers as character throughout melting. Add regression
cases for leading zeros, decimal-looking identifiers and scientific notation, for
both sample and taxon identifiers.

Location: `R/plot-methods.R`, `psmelt()`.

## 4. Sample merging destroys categorical metadata labels

**Priority: high for API design; confirmed behavior, partly intentional.**

Sources: [#608](https://github.com/joey711/phyloseq/issues/608),
[#1205](https://github.com/joey711/phyloseq/issues/1205),
[#1345](https://github.com/joey711/phyloseq/issues/1345).

`merge_samples()` coerces compatible metadata columns to numeric before
aggregation. In the local reproduction, factor values `soil` and `water` become
numeric values `1` and `2`, without a warning. Averaging conflicting factor codes
can also produce numbers with no categorical interpretation.

This is present in the implementation by design, so changing it requires a clear
compatibility policy. Preserve categorical values when constant within a group;
report conflicts explicitly and require a strategy for resolving them.

Related [#465](https://github.com/joey711/phyloseq/issues/465): `fun` being ignored
for OTU counts is explicitly documented; it should not be treated as an undisclosed
implementation bug. A new API could separate count aggregation from metadata
aggregation to remove the ambiguity.

Location: `R/merge-methods.R`, `merge_samples()` methods.

## 5. Single-column metadata and taxonomy break merging

**Priority: medium. Both confirmed locally; suitable early fixes.**

- [#1744](https://github.com/joey711/phyloseq/issues/1744): merging objects with
  one metadata column fails with `sam_data slot is empty`. The final data-frame
  subset drops dimensions. The issue discussion acknowledges the edge case.
- [#1359](https://github.com/joey711/phyloseq/issues/1359): `tax_glom()` with one
  taxonomic column fails with `dim(X) must have a positive length`. The failing
  subset is inside the `merge_taxa()` method for `taxonomyTable`.

Next change: preserve dimensions with `drop=FALSE`, with regression tests for
single rows, single columns, both OTU orientations and whole-object reconstruction.
The analogous one-sample OTU merge was also checked and works in this environment;
do not assume all visually similar subscripts are broken.

Location: `R/merge-methods.R`.

## 6. Filtering wrappers lose the caller's local variables

**Priority: medium. Both confirmed locally.**

Sources: [#1316](https://github.com/joey711/phyloseq/issues/1316),
[#487](https://github.com/joey711/phyloseq/issues/487).

Both `function(x, target) subset_taxa(x, Genus == target)` and the corresponding
`subset_samples()` wrapper fail with `object 'target' not found`. The subsetting
expression is evaluated through a wrapper without preserving the original
evaluation environment.

Next change: capture the expression and caller environment explicitly, while
maintaining column lookup precedence. Test local variables, nested wrappers and
name collisions. A base-R solution should be evaluated before adding dependencies.

Locations: `R/taxonomyTable-class.R`, `R/sampleData-class.R`.

## 7. Scaling and compatibility work

- [#1245](https://github.com/joey711/phyloseq/issues/1245),
  [#517](https://github.com/joey711/phyloseq/issues/517): `tax_glom()` repeatedly
  reconstructs objects while merging groups. Strong candidate for grouped
  aggregation, after establishing equivalence for taxonomy, reference sequences
  and trees. The reported large-dataset timings were not reproduced here.
- [#1041](https://github.com/joey711/phyloseq/issues/1041): `psmelt()` expands the
  complete taxon-by-sample matrix. The thread includes memory failures on large
  datasets. Benchmark allocation and consider an optional nonzero-only result;
  changing the default would change zero-row semantics.
- [#1794](https://github.com/joey711/phyloseq/issues/1794): `aes_string()` emits a
  deprecation warning under local ggplot2 4.0.3. Plot construction still succeeds.
  Migrate aesthetic mappings with coverage for names requiring quoting.
- [#1764](https://github.com/joey711/phyloseq/issues/1764): the reported sparse
  transpose crash was **not reproduced** with local Matrix 1.7-6. Test additional
  versions and load orders before changing S4 registration.

## Additional finding while checking rarefaction

Related workflow discussion:
[#1787](https://github.com/joey711/phyloseq/issues/1787).

The exact issue asks about repeated rarefaction; it does not report the following
failure. Independently, `rarefy_even_depth(..., rngseed=42)` fails with
`object '.Random.seed' not found` if no RNG state exists. The implementation reads
the old seed before initializing it. Reproduced after removing `.Random.seed` in
a disposable R process. Address absence/presence of RNG state, restoration on
errors and deterministic results in the rarefaction refactoring.

## Already covered by the local first change

[#1262](https://github.com/joey711/phyloseq/issues/1262), multi-argument
`transform_sample_counts()` callbacks, succeeds with the current local
transformation implementation. The fix is local and not yet published upstream.

## Suggested implementation order

1. Protect UniFrac results and stop silent losses in merging and `psmelt()`.
2. Fix the two single-column failures and filtering evaluation environments.
3. Define categorical metadata and incomplete-component policies with compatibility
   tests, then refactor merging around those policies.
4. Benchmark grouped aggregation and melting; migrate plotting compatibility.
5. Add normalization and contaminant workflows on top of the validated operations.

The original audit added only diagnostics. The 2026-11-08 changes address items
1 and 2 above; items 3–7 and the rarefaction RNG finding remain open.

## Reproducing the local audit

From the repository root, with development dependencies installed:

```sh
Rscript --vanilla .github/diagnostics/issue-audit.R
Rscript --vanilla .github/diagnostics/unifrac-audit.R
```

These are diagnostic scripts, not a passing regression suite: they catch errors
and report outputs so every case can be inspected. They intentionally exit
successfully even when a reported bug is present. Run them in fresh R processes;
the general diagnostic deliberately removes RNG state to inspect the rarefaction
failure. They load the local source without installing over another phyloseq.
