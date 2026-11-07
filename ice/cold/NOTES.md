# Cold adaptation in ice vs water (CMO mesocosms), 2026-10-24

Workstream: cold-adaptation genes and proteome indices in genomes and transcriptomes, plus HGT of cold genes. Population genetics and ice physics live elsewhere (ice/popgen, parent).

## Design and caveats
- **Paired contrasts:** RNA samples are paired within season × pool × day.
  - Ice vs water: 24 blocks (both seasons, both pools).
  - Ice TOP vs BOT: 13 blocks (2024-25 only).
- **Statistics are descriptive.** Each season has one pool per treatment, and samples within a pool are pseudoreplicates. "n/24" counts blocks, not independent replicates.
- **Expression units:** CPM of CDS read pairs. Gene IDs are BAKTA_EXTRA. KOs are the coordinate-corrected ones (data/genemap.pkl).
- **Pfam search:** hmmsearch `--cut_ga` of 17 cold-related Pfam-A families against all 13.3M proteins (Grex job 8120836, `grex/pf.domtbl`).
- **MAG ice/water preference:** taken from the parent's `ice/out/mag_icewater.tsv` (delta = median log2 ice/water over the 24 blocks, 396 MAGs).
  - Classes: ice-associated = delta ≥ 0.75 and q < 0.05 (44 MAGs); water-associated = delta ≤ −0.75 and q < 0.05 (24 MAGs); the remaining 328 are neutral.
  - My own simpler RNA+DNA version is in `mag_ice_assoc.tsv`; it is superseded.

## 1. Transcription of cold genes: ice vs water and TOP vs BOT
Columns: category; ice vs water (median log2 over blocks, blocks up); TOP vs BOT; within-MAG (category share of each MAG's own transcripts, ice vs water; median of per-MAG medians and % of MAGs up).

| category | ice vs water | TOP vs BOT | within-MAG |
|---|---|---|---|
| Ice-binding proteins | +2.10 (23/24) | +0.55 (9/13) | +1.96 (71% of 17) |
| Cold-shock proteins | +0.71 (20/24) | +0.24 (10/13) | +0.70 (84% of 167) |
| Cold DEAD-box helicases | +0.66 (20/24) | +0.19 (10/13) | +0.58 (77% of 167) |
| Oxidative stress | +0.42 (15/24) | +0.43 (10/13) | +0.05 (52% of 197) |
| Cold-acclimation translation/RNA | +0.39 (21/24) | +0.06 (10/13) | +0.24 (73% of 218) |
| Chaperones | +0.33 (14/24) | +0.05 (9/13) | +0.42 (71% of 236) |
| Compatible solute synthesis | +0.32 (20/24) | +0.34 (11/13) | +0.23 (67% of 158) |
| Compatible solute transport | +0.16 (14/24) | +0.25 (10/13) | −0.06 (48% of 89) |
| Fatty-acid desaturases | +0.13 (15/24) | +0.19 (12/13) | +0.09 (55% of 148) |
| EPS export | +0.10 (14/24) | +0.36 (11/13) | −0.09 (48% of 80) |
| DMSP catabolism | +0.05 (14/24) | +0.17 (10/13) | −0.10 (47% of 15) |
| Osmotic channels/K+ uptake | −0.03 (11/24) | −0.06 (6/13) | −0.10 (42% of 94) |
| Membrane FA elongation (fabF/A/B) | −0.04 (11/24) | −0.02 (5/13) | −0.06 (48% of 164) |
| Carotenoid synthesis | −0.09 (10/24) | +0.18 (10/13) | −0.06 (49% of 49) |
| PUFA synthase (pfa) | −0.57 (4/24) | −1.57 (2/13) | −0.35 (43% of 7) |

- **The cold-shock program is regulation, not community turnover.**
  - cspA, the cold DEAD-box helicases (deaD/csdA +0.64, rhlE +0.61, srmB +1.06), RNase R (+0.44, 22/24), NusA, IF2, trigger factor and PNPase rise in ice.
  - The same rise appears inside 73–84% of individual MAGs.
- **Compatible solutes.**
  - Ectoine synthesis rises in ice: ectB +0.71, ectC +0.82, ectD +1.04.
  - So do proline synthesis (proB +0.55, 22/24) and the opuC osmoprotectant transporter (+1.07).
  - Top ice > bottom ice for opuC/opuBD, betT, trehalose otsB/treS, and desaturases (12/13 blocks). Top ice is the colder, saltier-brine layer.
- **Oxidative stress.**
  - katE, ahpF and osmC rise in ice (+0.7 to +1.5) and at the top (+0.7 to +0.8).
  - Within MAGs the change is flat, so this is mostly carried by community composition.
- **PUFA synthase (pfa) falls in ice and especially at the top (−1.6).** Its few carriers (Shewanella lorbin_049 and comebin_1103_sub, Colwellia) are not top-ice taxa.
- **Ice thickness (2024-25 ice, n = 24 samples).**
  - Cold DEAD-box helicases ρ = +0.61 (p = 0.002), translation/RNA ρ = +0.58 (p = 0.003), cold-shock proteins ρ = +0.40 (p = 0.05).
  - IBPs show no trend (ρ = 0.00).
  - Thickness is confounded with time.

## 2. Proteome indices (amino-acid composition)
- **Expressed proteome is not measurably more cold-adapted in ice.**
  - Transcript-weighted composition per sample: ice − water differences are ≤0.2 sample-SD for every index.
  - Arg/(Arg+Lys) is −0.007 (lower in 18/24 blocks), GRAVY −0.011 (lower in 18/24), Gln+Thr −0.0015 (lower in 19/24).
  - The fixed-MAG-composition counterfactual gives the same picture, so neither turnover nor regulation shifts bulk composition.
- **Genome level (396 MAGs; Spearman after GC + family fixed effects, n = 292).**
  - Ice-associated MAGs have more polar uncharged residues (S/T/N/Q; ρ = +0.24, q < 0.001) and more Gln+Thr (ρ = +0.18, q = 0.009).
  - They have less Gly (ρ = −0.27, q < 0.001), which is opposite to the textbook "more Gly" psychrophile signature.
  - Arg/(Arg+Lys) is weakly lower (ρ = −0.11, p = 0.05, q = 0.2).
  - Genomes are larger (ρ = +0.19, q = 0.009) and chaperone gene density is lower (ρ = −0.16, q = 0.02).
  - Raw ice-associated MAG medians: Arg/(Arg+Lys) 0.43 vs 0.54 neutral and 0.43 water, but that difference mainly reflects GC and taxonomy (Flavobacteriia, Colwellia).
- **Within-genus pairs (10 genera with an ice-leaning and a water-leaning MAG differing by ≥ 1 log2):** no index differs consistently. Arg/(Arg+Lys) is lower in the ice MAG in 7/10 pairs (median −0.012, p = 0.16).

## 3. Ice-binding proteins (DUF3494 and related; 1,076 genes)
- **The strongest ice marker in the dataset.** IBPs are 4.3× higher in ice than water (23/24 blocks) and 1.5× higher at the top than the bottom (9/13). The rise holds within MAGs too.
- **Genome content.**
  - IBP gene density is the clearest genomic correlate of ice preference: ρ = +0.34 after GC + family, q < 0.001.
  - 38.6% of ice-associated MAGs carry ≥ 1 IBP, against 11.0% of neutral and 16.7% of water-associated MAGs (odds ratio 4.9, Fisher p = 1.8e-5).
  - Examples: CAJQLB01 semibin_021 (6 IBPs), Shewanella vesiculosa / psychromarinicola (2–3), Polaribacter vamb_3778 (3), Croceibacter (3), Cognaticolwellia metabat_8244 (2), Nonlabens comebin_3710 (1).
- **Who transcribes IBPs in ice.** Domains by tiara:
  - Binned bacteria 43%.
  - Eukaryotes 31% (141 genes): mostly centric diatoms (Thalassiosira/Chaetoceros, MarFERReT best hits of co-located proteins; 28 CPM in ice vs 6 in water), plus a choanoflagellate-like contig (Acanthoeca hits; 14 CPM ice vs 0.7 water).
  - Unbinned prokaryotes 15%.
  - The top single genes are Cognaticolwellia metabat_8244 (JPFHOD_000395048/050, 10–13 CPM in ice, ≈8× water) and a Thalassiosira-contig IBP (JPFHOD_011124790, 12 CPM, 8×).
- **Marinomonas primoryensis lorbin_161 (ice-enriched, delta +1.25, 21/24 blocks) carries an MpIBP/MpAFP homolog.**
  - It is on contig_0024173. Diamond hits against the UniProt MpAFP fragments: A1YIY2 94.5% over 438 aa; A1YIY3 97% / 89% / 64% over 66–146 aa pieces.
  - The locus is broken into 8 short ORFs (JPFHOD_002267650–657, ~1.8 kaa) next to a "retention module" protein. The native protein is ~1.5 MDa with ~120 RII repeats; here they are collapsed and frameshifted, so this is an assembly artefact, not a short gene.
  - Transcripts are ≈0.3 CPM in both ice and water: present but not induced.
  - An intact large RTX-like adhesin set (1.4–2.2 kaa ORFs, S-layer/Ig-like, IS4/IS5 adjacent) sits in syntenic flanks in M. arenicola lorbin_050. It is not expressed.

## 4. HGT of cold genes
- **IBPs are rarely mobile.**
  - 13.8% of IBP genes are mobile-associated (MGE-borne or IS ≤ 10 kb); 4.2% are MGE-borne (3 strict plasmids, 37 phage-contig calls).
  - Those genes carry 6.2% of ice IBP transcription (0.7% for MGE-borne).
  - In binned bacteria, 20% sit within 10 kb of an IS element.
- **Sharing funnel.**
  - 153 IBP pairs share ≥ 90% aa / 80% coverage across different contigs: 75 same genus, 65 cross-genus, 8 cross-phylum, 2 cross-domain, 3 unassigned.
  - Flank test (blastn of ±20 kb windows): 7 transfers have a shared block in unrelated backgrounds and 1 more is HGT-like (block 1.8–6.8 kb, 91–100% aa).
  - Best supported is Nonlabens lorbin_115 (Bacteroidota) ↔ GCA-001735895 lorbin_446_sub (Pseudomonadota): 99% aa, 6.8 kb block, no flank similarity.
  - Others: Polaribacter metabat_7088 ↔ Labilibaculum (2.4 kb); Polaribacter ↔ Cyclonatronum (Balneolota, 4.3 kb); Paraglaciecola comebin_757_sub ↔ Colwellia (2.5–3.4 kb); Massilia ↔ Actimicrobium; Runella ↔ Fulvitalea; Flavobacterium ↔ Aquirufa.
  - 9 more "contained" pairs (the shorter window is fully identical) are probably the same genome split across contigs.
  - None of the transferred IBP copies is expressed above 0.1 CPM.
- **Eukaryotic IBPs show ancient, not recent, transfer.** 115 of 141 have a bacterial IBP hit, at median 54% aa (max 97%, a single unresolved pair). This fits the known bacterial origin of diatom IBPs, without recent exchange in the pools.

## Files (all under ice/cold/)
- `scripts/`:
  - aacount.py, bingc.py, cold.sh: Grex amino-acid counts + Pfam hmmsearch.
  - cats.py: category definitions. icecommon.py: pairing.
  - cold_genes.py, cold_within.py, indices.py, mag_test.py, ibp.py, ibp_hgt.py, euk_ibp_tax.py, mpibp.sh, ibp_extract.sh, build_json.py.
- `grex/`: pf.domtbl, counts.npy + ids.txt (per-protein AA counts, ACDEFGHIKLMNPQRSTVWY), bin_gc.tsv, ibp_* extracts, mpibp.m8, MarFERReT/taxdump helpers.
- **Tables:**
  - cold_gene_members.tsv, cold_inventory_mag.tsv, cold_cpm.tsv, cold_paired.tsv, cold_ko_paired.tsv, cold_within_mag.tsv / cold_within_summary.tsv
  - mag_indices.tsv, sample_indices.tsv, index_paired.tsv, mag_cold_table.tsv, mag_test.tsv, genus_pairs.tsv
  - ibp_genes.tsv, ibp_taxa.tsv, ibp_domain.tsv, ibp_context.tsv, ibp_pairs.tsv, ibp_euk_nearest.tsv, ibp_euk_taxa.tsv, pf_hits.tsv, pf_dom.tsv
- `cold.json`: site data (see the README section in the parent's site build).
