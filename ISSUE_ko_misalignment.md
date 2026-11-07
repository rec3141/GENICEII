# KO / eggNOG / CAZy assignments are joined to the wrong genes (mag_analysis → illumina_rna)

Found 2026-10-23 while exploring the CMO metatranscriptomes. Not fixed in any pipeline code yet.

## Symptom
- In the marine run, `K00496` (alkane 1-monooxygenase) landed on genes Bakta calls RecA, tRNA pseudouridine synthase and HisF, among others.
- Genes Bakta calls "Alkane 1-monooxygenase" carried unrelated KOs.
- The pipeline's `ko_tpm`, `cog_tpm` and `cazy_tpm` tables are therefore wrong, and so are the Function tab of the public RNA viz (microscape.app/marine-coassembly-flye-rna).

## Root cause
- In `mag_analysis/main.nf`, every protein-level tool consumes **BAKTA_BASIC** outputs (`ch_proteins`, `ch_gff`):
  - KOFAMSCAN, EMAPPER, DBCAN, per-gene Kaiju, IslandPath, MacSyFinder, DefenseFinder, CALCULATE_GENE_DEPTHS, and MAP_TO_BINS.
- **BAKTA_EXTRA** is a second, independent Bakta run whose output only replaces `ch_annotation_tsv`. It numbers its locus tags (`JPFHOD_…`) independently, so the IDs agree near the start of the assembly and then drift.
- `illumina_rna/main.nf` `referencesFromMag()` picks `annotation/bakta/extra/annotation.gff3` first. Its genes therefore use extra IDs, while the `annotations` (merged_annotations.tsv) and `kofam` inputs it joins to are keyed by basic IDs.
- `illumina_rna/bin/augment_viz_functions.py` reads `bakta/basic/annotation.tsv`, which is consistent with the annotations but not with the GFF used for counting.

## Evidence
- Remapping by identical (contig, start, end) coordinates fixes the joins. 13,323,490 of 13,324,269 CDS match (99.99%), and K00496 then lands on "Alkane 1-monooxygenase" genes.
- Remap job: Grex `~/scratch/project_rnaseq/explore/koremap.sh`; output `data/ko_remap.tsv.gz` here.

## Fix options
1. **illumina_rna, minimal.** Prefer `bakta/basic/annotation.gff3` in `referencesFromMag()` (swap the two paths).
   - Gene IDs then match every mag_analysis protein annotation.
   - Cost: product names come from the light-db run, with more "hypothetical protein".
   - Extra-run products could still be attached by coordinates in `summarize_expression.py`.
2. **illumina_rna, robust.** Keep the extra GFF, but key annotations by coordinates.
   - Join kofam/eggNOG/dbCAN to genes through `bakta/basic/annotation.tsv` (ID → contig:start:stop), instead of by ID.
   - Also add a sanity check: on a random sample of genes, warn if Bakta product words and KO definitions rarely agree.
3. **mag_analysis, root.** When `--bakta_extra` is set, feed BAKTA_EXTRA proteins/GFF to all downstream protein tools.
   - Alternatively, make BAKTA_EXTRA reuse the basic run's gene calls (`--regions`/`--proteins`) so locus tags are identical.
   - This makes the per-MAG KEGG modules, viz and defense/MGE tables internally consistent too. Today MAP_TO_BINS uses the basic GFF, which is consistent with the annotations, but the viz mixes the extra TSV with basic-keyed annotations.

## After fixing
- Rerun `rna-marine-grex.sh final` (summarize stage only) and redeploy the RNA viz.
- Check other mag_analysis runs that used `--bakta_extra` (chesterfield, freshwater, per-genome mag_* runs) for the same split.
