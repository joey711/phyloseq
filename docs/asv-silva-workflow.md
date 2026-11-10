# Targeted 16S ASV phylogeny with versioned SILVA references

Implemented 2026-11-08. The package remains named `phyloseq`. This workflow
selects observed ASVs, filters bacterial/archaeal references and estimates a
small joint maximum-likelihood tree. Full reference genes provide context; they
do not reconstruct unobserved ASV bases or establish species identity.

## Inputs and dependencies

Use a phyloseq with `refseq`, a TreeSummarizedExperiment with `referenceSeq`, or
a named character/DNAStringSet of observed sequences. Sequence IDs must match
all container feature IDs exactly. IDs such as `01` remain characters; duplicate
DNA with different feature IDs remains distinct. ASVs must be ungapped A/C/G/T;
case is normalized and U becomes T. References may contain IUPAC ambiguities.

```r
pkgload::load_all(".") # development checkout, fresh R session
# Optional dependencies used only by infer_asv_tree():
# BiocManager::install(c("DECIPHER", "phangorn"), update=FALSE)

observed <- asv_sequences(ps)
selected <- select_asvs(ps, ids=c("ASV_1", "ASV_2"))
# Alternatively: select_asvs(ps, sequences=c("ACGT...actual DNA..."))
```

Replace example IDs with IDs in your own object. The sequence selector uses
literal complete observed sequences, without dots or wildcards. Set
`both_strands=TRUE` to include reverse complements. Selection preserves samples,
counts of retained features, annotations and existing assay provenance. Removing
features can make an old CLR/proportion assay stale: retransform from its source
before using that assay. `prune=FALSE` returns a named logical mask.

## Explicit release and download

```r
raw <- silva_reference(
  version="144", cache_dir="~/reference-databases/silva", download=TRUE
)
# Later, entirely offline with the same persistent cache:
raw <- silva_reference("144", cache_dir="~/reference-databases/silva")
```

Official SSU Ref NR99 downloads support releases `144` and `138.2`, using explicit
release URLs, FASTA and quality-file MD5 sidecars. The version is required;
there is no changing implicit `latest`. `download=FALSE` never accesses the
network. The temporary-directory default cache is session-lived: choose a
persistent cache for real use. License/citation files accompany explicit
downloads; retain attribution with exported analyses. The release 144 unaligned
FASTA alone is about 309 MB compressed, before the quality file and filtered
outputs. See the [official export listing](https://www.arb-silva.de/current-release/Exports).

For existing or other historical exports, supply both matching files:

```r
raw <- silva_reference("138.2", fasta="/data/reference.fasta.gz",
                       quality="/data/reference.quality.gz")
```

Local files receive computed checksums; `publisher_verified=FALSE` distinguishes
them from cache files checked against publisher sidecars. Checksums detect
changed contents; they are not an authenticity guarantee for a user-edited cache.
The selected release must match the files supplied by the user.

## Auditable reference QC

```r
filtered <- filter_silva_reference(
  raw, output_dir="silva144-16S-profile1",
  domains=c("Bacteria", "Archaea"),
  min_length=c(Bacteria=1200, Archaea=900), max_length=3000,
  max_ambiguous=0.02, max_homopolymer=0.02, max_vector=0.02,
  min_sequence_quality=95, min_alignment_quality=75,
  missing_quality="exclude"
)
filtered$summary
qc <- read.delim(filtered$qc_report)
table(qc$Keep)
# Reuse the saved descriptor (paths must still exist):
filtered <- readRDS("silva144-16S-profile1/reference.rds")
```

Thresholds above are editable defaults, not a universally optimal quality
profile. Eukaryotic 18S is excluded. Domain minima follow the SILVA Ref
distinction, but length and the `_trunc` filename do not prove complete genes.
The [SILVA release notes](https://www.arb-silva.de/documentation/release-144)
explain that export and the underlying reference filters.

QC requires sequence/alignment scores and ambiguity/homopolymer/vector metadata,
joined by **accession.start.stop**, including multiple rRNA regions on one
accession. Scores are 0–100; function percentage thresholds are fractions.
Computed and reported ambiguity must both pass. The homopolymer criterion uses
SILVA's percentage metric, not maximum run length. No per-read Phred scores are
invented from a FASTA. Missing required values are excluded with a reason, or
stop the run with `missing_quality="error"`. Invalid/unavailable requested fields
are errors. Optional `min_pintail` works only when the historical export supplies
that field; it errors for release 144, which discontinued Pintail. With the
default `min_pintail=NULL`, no additional chimera screening is claimed.

The source remains unchanged. A fresh output profile contains a filtered FASTA,
gzipped per-sequence QC table (including all exclusion reasons) and an RDS
descriptor with checksums, thresholds and counts. Failure counts can overlap.
An all-excluded result is inspectable but cannot support inference. Existing
output files are never overwritten. FASTA and report are streamed in chunks;
quality metadata are loaded in memory, so memory is not constant in database size.

## Candidate retrieval and targeted inference

```r
hits <- match_asvs_to_reference(selected, filtered,
                                max_mismatches=0, both_strands=TRUE,
                                max_hits=100)
hits$summary
head(hits$matches)

# Proceed after reviewing no-hit, ambiguous-orientation and truncation reports.
# An explicit subset is useful to choose a defensible local reference context:
# reference_ids <- c(...) # actual IDs from hits$matches$ReferenceID
result <- infer_asv_tree(hits, model="GTR", gamma_categories=4,
                         bootstrap=100, seed=42, threads=1)
plot(result$tree, cex=0.7)
result$tip_data
result$provenance
```

Matching uses Biostrings' full-query **ungapped** search, exact by default, with
an explicit integer mismatch allowance. Identity is matching sites divided by
query length; QueryCoverage is 1. A short query need not cover the full reference.
This first matcher is appropriate for targeted selections, not a benchmarked
replacement for indexed VSEARCH/BLAST on thousands of ASVs. Indels can cause
true homologs to have no hit: use a future gapped retrieval backend rather than
claiming such queries lack relatives. Manually curated named reference sequences
can also be passed, explicitly bypassing SILVA QC.

Multiple references and reverse complements remain explicit. Candidate caps are
reported as `Truncated`; no-hit queries remain in the summary. Inference stops on
no-hit queries, unresolved orientation, fewer than three references, insufficient
pairwise observed-site overlap or more than `max_sequences=300` total sequences.
For truncation, increase `max_hits` or deliberately pass `reference_ids`; every
query must retain a candidate. This is a user-chosen subset, not recovery of the
omitted uncertainty. Ambiguous strand requires a named `query_strands` vector of
`+`/`-` for every query. No species confidence or placement probabilities are
derived from these matching scores.

DECIPHER aligns the actual fragments and reference genes together; phangorn
fits JC or GTR, optionally with discrete gamma and nonparametric ML bootstrap.
The seed is explicit and the caller's RNG state is restored. Alignment must
preserve each ungapped sequence; branch lengths must be finite and nonnegative.
The fitted model, aligner versions, IDs, sequence fingerprints and matching/QC
settings are retained. A zero bootstrap means **support was not assessed**.
The result does not supply phylogenetic-placement uncertainty. Biological
accuracy and sensitivity to fragment length/reference selection still need
validation on appropriate biological datasets.

## Output and analysis boundaries

- `tree`: unrooted context tree, labelled `REF::<ID>` and `ASV::<ID>`.
- `tip_data`: explicit mapping to reference IDs and observed feature IDs.
- `asv_tree`: pruned tree using original ASV IDs, or NULL for one ASV. Full-context
  bootstrap labels are removed because they are not ASV-only clade support.
- `alignment`, `fit`, `bootstrap_trees`, `query_sequences` and `provenance`:
  inspectable results; observed sequences and input counts remain unchanged.

The method is a **joint ML tree**, not placement on a fixed reference backbone.
EPA-ng/profile alignment, native jplace uncertainty and a scalable indexed search
are future integrations. Full-length reference context cannot increase the
number of sites observed in an amplicon. Disjoint target regions need a more
appropriate placement workflow.

No tree is automatically attached to phyloseq/TSE. The context tree contains
reference tips that are not counted features. Before downstream UniFrac, make
an explicit scientifically justified rooting decision on a suitable ASV tree
and verify tip/count agreement; this workflow does not invent an outgroup.
The legacy UniFrac wrapper's behavior on unrooted inputs is unchanged.

## Reproducible offline demonstration and validation

```sh
Rscript --vanilla .github/diagnostics/asv-silva-demo.R /tmp/asv-silva-demo
Rscript --vanilla -e 'testthat::test_local(".", filter="asv-silva|phylogeny-integrity", reporter="summary", stop_on_failure=TRUE)'
```

The demo creates synthetic full genes, ASV fragments and quality metadata,
executes the local-file workflow and saves an alignment, candidate tables, QC,
RDS result and annotated tree PNG. It downloads nothing and is not a biological
benchmark. Tests cover selection/containers, ID collisions, duplicate DNA,
reverse complements, QC boundaries, missing/invalid fields, cache checksums,
candidate ambiguity/truncation, ML/bootstrap integrity and RNG restoration.
The complete release-144 download and real-export integration check now pass:
905,628 input references, 267,459 retained under the default QC profile, full
filtered-database search, and a targeted 13-tip ML tree. See the
[live validation report](validation/silva144-live.md). Release 138.2 has header,
checksum and offline schema coverage; its complete download has not been tested.

For an explicit, larger integration check with official release-144 files:

```sh
Rscript --vanilla .github/diagnostics/validate-silva-live.R /path/to/cache /path/to/new-output --download
# Reuse the same downloaded cache and saved QC result offline:
Rscript --vanilla .github/diagnostics/validate-silva-live.R /path/to/cache /path/to/new-output
python3 .github/diagnostics/verify-silva-files.py /path/to/cache /path/to/new-output
```

Only the `--download` flag permits network access. This check applies QC to the
entire database, searches two known-origin 250-nt reference-derived fragments,
then fits a targeted joint tree using up to ten candidates per fragment. The
context choice is a deterministic integration fixture, not a recommended
biological reference-selection strategy. It saves timings, checksums, QC counts,
candidate tables, alignment, fitted result and session information outside the
repository. This validates integration with real exports, not accuracy on
independently observed ASVs or placement uncertainty.

The two fragments use overlapping coordinates on one parent gene. Equal raw
offsets on different partial genes would not guarantee a homologous shared
region. The independent Python check streams the original/filtered FASTA headers
and QC rows, comparing counts and exclusion consistency without calling R.

Domain selection uses SILVA's taxonomy. Chloroplast/mitochondrial entries can
occur within bacterial taxonomy; this profile does not implement a separate
organelle-exclusion criterion. Domain/QC filtering is not general contaminant
detection or proof that a sequence originates from a free-living prokaryote.
