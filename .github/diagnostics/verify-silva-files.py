"""Count real FASTA/QC files independently of R. No network or dependencies.

python3 .github/diagnostics/verify-silva-files.py CACHE OUTPUT
"""
import collections
import csv
import gzip
import json
import sys
from pathlib import Path

if len(sys.argv) != 3:
    raise SystemExit("Usage: verify-silva-files.py CACHE OUTPUT")
cache, root = (Path(value).expanduser() for value in sys.argv[1:])


def check(condition, message):
    if not condition:
        raise ValueError(message)


source = cache / 'SILVA_144_SSURef_NR99/SILVA_144_SSURef_NR99_tax_silva_trunc.fasta.gz'
domains = collections.Counter()
with gzip.open(source, 'rt') as stream:
    for line in stream:
        if line.startswith('>'):
            taxonomy = line.strip().split(None, 1)[1]
            domains[taxonomy.split(';', 1)[0]] += 1
qc_domains = collections.Counter()
retained_domains = collections.Counter()
reasons = collections.Counter()
records = kept = 0
with gzip.open(root / 'qc/reference-qc.tsv.gz', 'rt') as stream:
    for row in csv.DictReader(stream, delimiter='\t'):
        records += 1
        qc_domains[row['Domain']] += 1
        check(row['Keep'] in ('TRUE', 'FALSE'), 'Invalid QC keep flag')
        if row['Keep'] == 'TRUE':
            check(not row['Reason'], 'Retained record has exclusion reasons')
            retained_domains[row['Domain']] += 1
            kept += 1
        else:
            check(bool(row['Reason']), 'Excluded record has no reason')
            reasons.update(row['Reason'].split(';'))
with (root / 'qc/reference-16S.fasta').open() as stream:
    fasta_records = sum(line.startswith('>') for line in stream)
check(domains == qc_domains, 'Source/QC domain counts differ')
check(sum(domains.values()) == records, 'Source/QC total counts differ')
check(fasta_records == kept, 'Filtered FASTA/QC retained counts differ')
result = dict(source_records=sum(domains.values()), source_domains=dict(domains),
              qc_records=records, retained=kept, excluded=records-kept,
              retained_domains=dict(retained_domains), reasons=dict(reasons),
              filtered_fasta_records=fasta_records)
(root / 'independent-counts.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
