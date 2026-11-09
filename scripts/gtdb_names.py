"""Give every refined DAS Tool bin a GTDB row under its plain name (rerunnable).

data/gtdbtk_taxonomy.tsv is GTDB-Tk r226 run on all binner outputs plus the DAS Tool set. A DAS Tool bin that is
an unchanged binner bin appears twice (`dastool-lorbin_184` and `lorbin_184`, identical rows). A bin that DAS Tool
refined (`*_sub`) exists only as `dastool-comebin_4318_sub`. The site's scripts look GTDB up by the plain name, so
those bins fell back to Kaiju or "unclassified".

Rule: for every `dastool-X` row with no `X` row, append a copy named `X`. Rows that exist are never changed, so a
second run adds nothing. Run it again whenever data/gtdbtk_taxonomy.tsv is copied from Grex.
"""
import sys
import pandas as pd

F = '/data/dev/rnaseq-exploration/data/gtdbtk_taxonomy.tsv'
gt = pd.read_csv(F, sep='\t', dtype=str, keep_default_na=False)
names = set(gt.user_genome)
ds = gt[gt.user_genome.str.startswith('dastool-')]
plain = ds.user_genome.str[len('dastool-'):]

# check: where both names exist the rows must agree, otherwise the alias rule is not safe
both = ds[plain.isin(names)].assign(p=plain[plain.isin(names)]).merge(gt, left_on='p', right_on='user_genome', suffixes=('', '_p'))
bad = both[both.classification != both.classification_p]
if len(bad):
    sys.exit(f'{len(bad)} bins have different GTDB calls under the dastool- and plain names; not adding aliases')

add = ds[~plain.isin(names)].copy()
add['user_genome'] = plain[~plain.isin(names)]
if len(add):
    pd.concat([gt, add]).to_csv(F, sep='\t', index=False)
print(f'{len(ds)} dastool- rows, {len(both)} already had a plain-name row (all identical), added {len(add)} plain-name rows '
      f'({add.user_genome.str.endswith("_sub").sum()} *_sub)')
