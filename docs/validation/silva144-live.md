# SILVA 144 real-export integration check — 2026-11-08

**Passed:** complete official download, publisher MD5 checks, whole-database QC,
full filtered-database candidate search, targeted joint ML inference, and an
independent streamed comparison of source/report/filtered-file counts.

This is an integration check using reference-derived fragments. It does not
establish accuracy on independently observed ASVs or validate backbone placement.
Runtime package code and scientific QC/overlap thresholds were unchanged during
this retry; the additions are opt-in diagnostics and this report.

## Source and profile

Dataset: explicit release **144, SSU Ref NR99**, unaligned gene-trimmed FASTA
and its matching quality export. The original downloaded files remain unchanged.

| File | Publisher-verified MD5 |
| --- | --- |
| `SILVA_144_SSURef_NR99_tax_silva_trunc.fasta.gz` | `ccead7bba52075f5d8ba2bfab5adc4c1` |
| `SILVA_144_SSURef_Nr99.quality.gz` | `840ac0bdc1efa1264035cc8576e020e9` |

The files and sidecars came from the explicit
[release-144 exports](https://www.arb-silva.de/archive/release_144/Exports),
including its `quality` directory. License/citation files were downloaded with
the cache. Cache checks and a subsequent offline run both passed.

The unchanged default profile selects Bacteria/Archaea, minimum lengths
1,200/900 nt, maximum 3,000 nt, maximum ambiguity/homopolymer/vector fractions
0.02, minimum sequence quality 95 and alignment quality 75. Missing required
quality values are excluded. No Pintail or additional chimera screening was
requested. These defaults are explicit choices, not universally optimal QC.

## Counts independently verified from actual files

| Domain | Original FASTA/QC rows | Retained FASTA/QC rows |
| --- | ---: | ---: |
| Bacteria | 626,354 | 256,126 |
| Archaea | 32,178 | 11,333 |
| Eukaryota | 247,096 | 0 |
| Total | **905,628** | **267,459** |

The profile excluded **638,169** records. Domain counts above were counted from
the downloaded FASTA taxonomy headers and independently compared with every QC
row. They are file-derived observations, not copied website summary categories.
The retained FASTA contains exactly the reported 267,459 sequence records.
Each excluded row has a reason; retained rows have none.

| Exclusion reason | Records failing criterion |
| --- | ---: |
| Domain | 247,096 |
| Sequence quality | 582,128 |
| Alignment quality | 60,742 |
| Length | 677 |
| Ambiguity | 1 |
| Missing quality / homopolymer / vector | 0 |

Reasons overlap and must not be summed as mutually exclusive categories. The
sequence-quality threshold is restrictive on these exports; inspect thresholds
and reports rather than treating retention count as a quality guarantee.

SILVA can classify chloroplast/mitochondrial entries within bacterial taxonomy.
This profile does not independently exclude organelles. It excludes eukaryotic
18S by domain but is not a general contaminant-removal or host-origin assay.

## Candidate search and tree

Two unambiguous 250-nt fragments were derived from one passing reference
(`HG374866.9.1343`), using positions 300–549 and 325–574. They share known parent
coordinates. This is a positive integration fixture, not independently observed
ASVs from a microbiome sample.

The exact, both-strand search scanned all 267,459 retained references:

| Query | Candidate references | Returned | Truncated | Strand |
| --- | ---: | ---: | --- | --- |
| `derived_01` | 85 | 85 | FALSE | + |
| `derived_02` | 86 | 86 | FALSE | + |

The fixture explicitly selects up to ten returned references per query, producing
11 unique reference genes plus two query fragments. This deterministic choice
tests the software path; it is not a biological diversity/outgroup-selection
strategy. DECIPHER/phangorn produced a **13-tip, 1,530-site joint GTR ML tree**
with four gamma categories, seed 42 and finite nonnegative branch lengths. The
fragment-only tree preserves both query IDs; returned query sequences are unchanged.
Bootstrap was zero: **support was not assessed**.

An initial fixture using equal raw offsets on two different genes stopped at the
overlap check. The fixture was corrected to shared parent coordinates. Equal
raw offsets on partial/different genes do not guarantee homologous overlap. The
production overlap guard was retained; no failed tree was accepted by lowering
its threshold.

## Local timings and validation scope

On the successful offline reuse run (R 4.4.3, macOS arm64):

| Stage | Elapsed seconds |
| --- | ---: |
| Verify cached source files | 0.975 |
| Search two queries across the filtered database | 29.884 |
| Targeted alignment and joint ML inference | 23.491 |

These are single-run diagnostic timings with downloaded cache and existing QC
outputs. They omit the initial download/full-QC time, which was not retained
when the initial fixture stopped. They are not a before/after benchmark or a
general performance guarantee. Backends: Biostrings 2.74.1, DECIPHER 3.2.0,
phangorn 2.12.1. Full package tests/check previously passed 1,310 expectations
without test warnings/skips; runtime package source was not changed in this retry.

Remaining validation includes interrupted downloads, a complete historical
release, real independently measured ASVs, biological accuracy, bootstrap
behavior on biological datasets and reference-backbone placement uncertainty.

## Reproduction and artifacts

From the repository root, choose a cache and output directory:

```sh
Rscript --vanilla .github/diagnostics/validate-silva-live.R /path/to/cache /path/to/output --download
Rscript --vanilla .github/diagnostics/validate-silva-live.R /path/to/cache /path/to/output
python3 .github/diagnostics/verify-silva-files.py /path/to/cache /path/to/output
```

Normal CI never downloads this database. The explicit flag permits downloads;
the second command reuses and verifies cached files/QC outputs offline. The
Python standard-library check reads original/filtered FASTA and QC tables without
calling the package.

This machine's files are outside Git:

- Cache: `/private/tmp/phyloseq-silva-live-cache`.
- Outputs: `/private/tmp/phyloseq-silva-live-validation`.
- `validation.rds`: reference fingerprints, counts, candidates, model settings,
  query/reference hashes, timings and session details.
- `independent-counts.json`: independent file-count and reason validation.
- `phylogeny.rds`, `alignment.fasta`, `candidates.tsv`, `query-summary.tsv` and
  `qc/reference-qc.tsv.gz`: inspectable outputs.

These are temporary paths and may be cleaned by the operating system. Use a
persistent directory for long-term work. No database or external binary is
bundled in the R package. No commit or push was made during this validation run.
