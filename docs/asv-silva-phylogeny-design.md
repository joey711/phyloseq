# ASVs, versioned SILVA references and phylogeny

Design and roadmap, updated 2026-11-08. A first targeted **16S bacterial/archaeal**
workflow is implemented: `asv_sequences()`, `select_asvs()`, `silva_reference()`,
`filter_silva_reference()`, `match_asvs_to_reference()` and `infer_asv_tree()`.
See the [usage guide](asv-silva-workflow.md) for the actual API and limitations.
The inference backend is a joint DECIPHER/phangorn ML tree; reference-backbone
placement, gapped candidate retrieval and additional chimera screening remain
future work. Validation includes offline fixtures, a synthetic demo and a complete
official release-144 download/QC/search with targeted ML inference. The live test
uses reference-derived fragments, not independently observed ASVs or a biological
accuracy benchmark. See the [validation report](validation/silva144-live.md).

## Recommended direction

Offer a workflow callable from phyloseqV2 that selects ASVs, obtains an explicit
SILVA release, filters reference quality, identifies suitable references, and
places ASV fragments in a reference phylogeny. Keep phyloseq and
TreeSummarizedExperiment as user-facing containers. Full reference sequences
provide phylogenetic context; an ASV match does not reconstruct its unobserved
full-length gene or establish a species identity.

The reproduced polytomy bug in `fastUniFrac()` and partial-reference-sequence
merge losses are now addressed, with independent branch-incidence and merge
regression tests. `merge_phyloseq()` stops by default on incomplete refseq;
explicit `refseq_policy="drop"` discards that component with a warning. Other
legacy components retain their existing reconciliation behavior. See the
[upstream issue audit](../.github/ISSUE_AUDIT.md), especially sections 1 and 2.

## Current SILVA findings that affect the design

- The current SSU release is 144, released 2026-08-26. Its Ref NR 99 contains
  905,628 entries. LSU was not released as 144; version availability must be
  checked separately for each subunit, rather than constructing unverified URLs.
- Release 144 discontinued Pintail values. A requested Pintail threshold must
  fail clearly for this release, rather than being ignored or treated as passed.
- SSU Ref already applies length and alignment filters, but inclusion is not a
  guarantee that every reference is complete or suitable for every ASV region.
  Minimum lengths differ by domain. Additional completeness and target-region
  coverage checks are needed for a strict full-length analysis.
- Separate quality and full-metadata exports are available, with checksums.
  For SSU Ref NR99 release 144, the site lists a roughly 309 MB compressed
  unaligned FASTA, a 15 MB quality export and a 2 GB aligned FASTA. Stream/index
  the source instead of requiring the entire alignment in memory.
- The `_trunc` export suffix means unaligned nucleotides were removed. It does
  not mean an amplicon-specific reference subset. Do not equate filename suffixes
  with biological completeness.

Sources: [release 144](https://www.arb-silva.de/documentation/release-144),
[exports](https://www.arb-silva.de/current-release/Exports),
[quality exports](https://www.arb-silva.de/current-release/Exports/quality),
[field definitions](https://www.arb-silva.de/fileadmin/silva_databases/current/Fields_description/SILVA_description_of_fields.html),
[export README](https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/README.txt),
[archive](https://www.arb-silva.de/archive).

## Longer-term interface and data flow

This table describes the broader placement architecture, not additional exports.
The first delivery combines targeted alignment and ML inference in
`infer_asv_tree()`; its separate outputs do not include placement distributions.

| Proposed operation | Responsibilities | Output |
| --- | --- | --- |
| `select_asvs()` | Select by feature ID or DNA sequence; validate refseq/row links, IDs, strand and duplicate-sequence policy | Selected container plus an explicit ID/sequence map |
| `silva_reference()` | Resolve an explicit version/subunit/dataset, download files and quality metadata, verify official checksums and cache them | Reference descriptor and immutable raw-file paths |
| `filter_reference_quality()` | Apply recorded, version-aware sequence and metadata criteria | Filtered reference descriptor plus per-reference exclusion report |
| `match_reference_sequences()` | Search ASVs against filtered references, report query coverage, identity, orientation and ambiguous matches | One-to-many ASV/reference candidate table |
| `build_reference_tree()` | Select a sufficiently diverse reference context, align references, fit a declared model and assess support | Reference alignment, model and tree |
| `place_asvs()` | Align fragments to reference coordinates and perform phylogenetic placement | Native placements, uncertainty diagnostics, display tree and optional ASV-only tree |

The retention policy for equal/near-equal reference matches must be explicit.
One short ASV can be compatible with several different full-length references.
Do not silently substitute the best-hit reference sequence for the ASV, merge
its counts with a reference, or claim that the unobserved bases were measured.

Select enough references to represent relevant neighboring diversity and a
defensible outgroup when needed. A per-ASV hit cap must be recorded along with
whether results were truncated. Choosing only one best hit per ASV can remove
the context needed to assess ambiguous placement.

## Reference quality

Reference FASTA does not contain per-read Phred scores. Distinguish observable
sequence properties from SILVA's sequence/alignment quality metadata.

The quality report should support:

- Gene/subunit, domain and organism-group selection; separate 16S and 18S analyses
  when their different alignment and reference requirements warrant it.
- Length limits appropriate to the target and domain, plus coverage of the
  query's region. A length threshold alone is not proof of a complete gene.
- Ambiguous-base fraction, homopolymer criteria and vector-contamination metrics,
  distinguishing percentages from homopolymer run length.
- Available SILVA `seq_quality_slv`, `align_quality_slv`, `align_ident_slv` and
  related alignment metrics. Verify units against the selected release schema.
- Optional additional chimera screening with its method and dependencies
  recorded. Removed/missing Pintail values cannot act as a chimera check, and
  reference catalog entries do not supply sample abundance for abundance-based
  chimera assumptions.
- Exact-duplicate handling with an accession map, plus optional type-strain or
  cultured-sequence preferences. Preferences must not silently override QC.

Thresholds need explicit, editable presets rather than an undocumented universal
"high quality" switch. Requested unavailable fields should cause an error;
missing values need an explicit exclude/error policy. Match metadata on full
reference identifiers (accession plus gene coordinates), not accession alone:
one genomic accession can contain several rRNA regions.

Retain the raw reference files unchanged. Store the filtered IDs, parameters,
reason codes and a content hash separately so that decisions can be audited and
the same downloaded release can support multiple quality profiles.

## Phylogenetic methods

For short amplicon ASVs, the recommended main path is phylogenetic placement
onto a reference backbone. EPA-ng accepts a reference alignment, reference tree,
query alignment and fitted model; its native jplace result retains placement
weights instead of only a single selected attachment.

An alternative joint maximum-likelihood tree containing references and actual
ASV fragments can be useful for small targeted analyses. Label that method
separately and assess the influence of fragment length and missing data. Full
references do not make short query sequences more informative than their
observed sites. Do not describe this as recovery of full-length ASVs.

Candidate backends:

- R/Bioconductor for container handling, sequence QC, indexed/chunked retrieval
  and reports: Biostrings, DECIPHER and optional cache/database support.
- VSEARCH or BLAST for scalable candidate retrieval, validating the denominator
  used for identity/coverage when matching short queries to long references.
- SINA or MAFFT fragment/profile alignment for a reference coordinate system.
- IQ-TREE or another explicit ML backend for a reference backbone/model.
- EPA-ng for placement; retain its native output and version/model details.

The reference alignment, mask, model, tree tips and query coordinates must agree.
A taxonomy hierarchy exported as Newick is not an inferred sequence phylogeny.
SILVA's curated guide tree is a useful resource, but its suitability, branch
lengths and model still need validation for a particular placement workflow.

Sources: [EPA-ng documentation](https://github.com/pierrebarbera/epa-ng),
[SILVA guide-tree and alignment notes](https://www.arb-silva.de/documentation/release-144).

## Integration and reproducibility

"Included in the package" should mean one R workflow with dependency checks,
parameter handling, progress and errors. Download reference databases into a
user-selected versioned cache rather than distributing a large database inside
the package. External binaries need explicit installation/version checks; an R
wrapper does not bundle those executables automatically.

Record the database release/subunit/dataset, source URLs and checksums, quality
profile and schema, query sequence hashes, ASV/reference mapping, aligner and
model settings, mask, seed, threads, program versions and tree/placement hashes.
Resolve a convenience "latest" request to a concrete release before starting;
saved results must always refer to that concrete release.

Keep distinct outputs:

1. A context/display tree containing reference and ASV tips, with their origin
   explicitly annotated.
2. An optional analysis tree whose tip IDs match exactly the selected ASVs in
   the count matrix. Do not introduce reference sequences as observed features.
3. Original placement distributions, unmatched/low-confidence ASV reports,
   quality exclusions and a reference-accession table.

Attaching a tree must not drop samples/taxa, alter counts, replace observed ASV
sequences with reference genes, or overwrite unrelated TSE assays. Explicit
rooting and finite nonnegative branch-length checks are needed before UniFrac;
do not silently add a random biological root or invent missing branch lengths.

The current SILVA data license is CC BY 4.0 with attribution requirements. Store
the license/citation corresponding to the chosen release; older archived
releases may have different historical terms. See the
[SILVA license page](https://www.arb-silva.de/silva-license-information).

## Feasibility on this computer

The development environment currently has Biostrings 2.74.1, DECIPHER 3.2.0,
phangorn 2.12.1, BiocFileCache 2.14.0 and IQ-TREE 2.3.5 available. VSEARCH,
BLAST, MAFFT, SINA and EPA-ng were not found on PATH during this inspection.
Availability is not a claim that their integration has been implemented or
tested. No external binaries were installed. The subsequent release-144
integration check downloaded its database into an external temporary cache,
as recorded in the validation report.

## Delivery priorities and acceptance criteria

The topology/refseq fixes, ASV selection, release resolution, streaming QC,
ungapped candidate mapping and targeted joint-tree backend are implemented.
Next priorities for this workflow are scalable gapped retrieval, a stable
reference-coordinate policy and explicit backbone placement with native jplace
uncertainty. The full release-144 download and integration path pass; interrupted
downloads, other complete releases and biological placement accuracy remain
additional validation work.

Acceptance coverage should include synthetic known placements, preserved
numeric-looking IDs, reverse complements, repeated DNA with different feature
IDs, one-to-many hits, no-hit queries, checksum failures, interrupted downloads,
unavailable releases/subunits/QC fields, missing quality values, disk-cache reuse,
container immutability, tree/count alignment and reproducibility. Use small
offline fixtures in normal CI; separate opt-in real-database and binary-backend
checks. Biological benchmarking is additional to software test success.

Other useful development priorities remain: explicit contaminant-control
integration, sparse/delayed assay processing and further compositional methods
with recorded zero/model assumptions. These should not delay the existing
data-integrity fixes.
