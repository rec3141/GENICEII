# HGT and diesel biodegradation in the CMO mesocosms (2026-10-23)

Both seasons were oiled with diesel: 2024-25 on DAY05 at 16:30, and 2025-26 (~136 L) on DAY01. There is one oiled (EXP) pool and one control (CTL) pool per season, so samples within a pool are pseudoreplicates and every statistic below is descriptive.

## Gene set (fixed after the KO misalignment)
The original genemap KO column came from eggNOG run on the Sep-24 *basic* Bakta annotation. Its JPFHOD locus tags are shifted relative to the *extra* annotation that the RNA reference uses. The gene set now rests on three sources, in this order:
1. Bakta extra products, which are self-consistent with the reference.
2. The coordinate-remapped KO table (data/ko_remap.tsv.gz).
3. CANT-HYD HMMs at the trusted cutoff.

Tiers:
- **confident:** 10,775 genes, supported by a product match and/or CANT-HYD.
- **putative:** 13,860 genes, supported by the remapped KO only. Many are broad orthology transfers onto hypotheticals, desaturases and VOC/cupin proteins.

Unless stated otherwise, numbers below use the confident tier. CANT-HYD contributes 337 CYP153 and 55 AlkB genes. Its AlkB cutoff is much stricter than the 577 product-named alkB.

## 1. Mobile context of the HC genes (confident, n = 10,775)
| flag | genes | % |
|---|---|---|
| geNomad plasmid, strict (score ≥ 0.9 with a hallmark or conjugation genes) | 139 | 1.3 |
| ...with conjugation genes | 109 | 1.0 |
| geNomad plasmid, any (score ≥ 0.7, includes false positives) | 1,462 | 13.6 |
| phage contig (≥ 2 hallmarks, score ≥ 0.95, CheckV host < viral) | 20 | 0.2 |
| provirus | 2 | <0.1 |
| IslandPath genomic island | 54 | 0.5 |
| integron (± 1 kb) | 20 | 0.2 |
| transposase/IS/integrase within 10 kb | 613 | 5.7 |
| **MGE-borne (any of the first set)** | 229 | 2.1 |
| **mobile-associated (MGE-borne or IS ≤ 10 kb)** | 734 | 6.8 |

- **Loose plasmid calls.** geNomad "plasmid_any" calls are unreliable on these long-read contigs. contig_0030642 (Cycloclasticus lorbin_184, dmp/tmo/meta-cleavage genes, plasmid score 0.94) is chromosomal: it is colinear with C. zancles and has a copy ratio of 1.1 (cyclo fork). Left in, it alone would make plasmids look like they carry 19% of the 2025-26 induction. It is excluded everywhere.
- **Viral calls.** geNomad's plain virus calls are very noisy (272k "virus" contigs). The strict phage set is tiny, and no phage carries an oil-induced HC gene.
- **GC.** Deviation from the contig GC does not separate mobile from chromosomal HC genes (mean ΔGC +0.8 for both).

## 2. Cross-taxon sharing (evidence for recent HGT)
Funnel:
- 6,217 pairs of HC proteins at ≥ 95% aa identity / ≥ 80% coverage on different contigs;
- 381 of those at ≥ 99% nt identity;
- 69 also in a different genus after GTDB/NCBI harmonisation, with uncorrelated DNA depth (ρ < 0.2 across 276 nanopore samples; same-MAG contig pairs have median ρ 0.60, random pairs −0.52);
- 26 flank-confirmed transfers: the identical block (0.6–7.5 kb) sits in ±20 kb flanks with no other alignment.

Of the 26, 9 are cross-class, 11 cross-family and 6 cross-genus. The rest of the 69 had flanks aligning at 90–97% throughout, i.e. related genomes rather than transfer.

Most "cross-family" sharing at 95% aa is artefactual or vertical:
- GTDB vs NCBI synonyms (Cycloclasticaceae vs Piscirickettsiaceae, Rhodobacteraceae vs Roseobacteraceae/Paracoccaceae).
- Redundant contigs of the same genome (100% nt with correlated depth).
- Congeneric strains (~94% nt with similar flanks).

Flank-confirmed HGT by category:
- **alkB (8 pairs):** one near-identical alkB block (99.2–100% nt) is shared by Imperialibacter (Bacteroidota), Hyphomonas, Tsuneonella (Erythrobacteraceae), Tardiphaga, Bdellovibrio and Pelagibacter/Fonsibacter contigs. The blocks are 1.6–4.9 kb, carry 2–6 genes, and sit in unrelated flanks.
- **Ortho-cleavage (cat/pca), 5 pairs:** Marimicrobium↔Marinobacter, UBA4592↔Pseudonocardia, Anderseniella↔Rhizobium, Paracoccus↔Pelagibacter, Psychromonas↔Vibrio.
- **RHO dioxygenase, 4 pairs:** e.g. Halioglobus↔Polaromonas, an 8-gene, 4.5 kb block.
- **Phenylacetate (paa), 4 pairs.**
- **BMM phenol hydroxylase operon:** dmpKLMNOP, 5 genes, 100% nt over 3.8 kb, between Hyphomicrobium (Alphaproteobacteria) and Rhodoferax (Betaproteobacteria).
- **alkG/alkT:** Sphingomonas↔Moraxella, 100%.
- Single pairs of extradiol, intradiol and meta-lower genes.

Most of the transferred copies are barely transcribed in the mesocosms: only 7 of 26 pairs reach ≥ 20 reads. The best expressed are:
- Halioglobus/Polaromonas RHO (70 and 147 reads, l2fc +1.2 to +1.8 oiled vs control, both seasons);
- the Pelagibacter/Fonsibacter alkB (about +3 in 2024-25);
- the Imperialibacter alkB (+2.7 to +3.0 in 2024-25).

## 3. Is the diesel response carried by mobile genes? Mostly not.
Share of the oil-induced HC transcript increase, Σ max(0, oiled − control) CPM over confident genes:

| | 2024-25 | 2025-26 |
|---|---|---|
| mobile-associated (MGE or IS ≤ 10 kb) | 9.2% | 8.3% |
| MGE-borne | 5.7% | 1.4% |
| strict plasmid | 5.4% | 0.5% |
| IS ≤ 10 kb | 8.4% | 8.3% |

- **Fold change.** Median l2fc (oiled vs control) of expressed genes (≥ 20 reads) is not different for mobile and chromosomal genes. 2024-25: 1.05 vs 0.70 (p = 0.50). 2025-26: 0.00 vs 0.33 (p = 0.20).
- **Who carries the induction.** Cycloclasticus accounts for 45.8% (2024-25) and 86.3% (2025-26) of the induced HC transcripts. Only 2.5% and 6.4% of that comes from mobile-associated genes: its 260 confident HC genes include 0 on plasmids, 1 on an island and 14 within 10 kb of an IS element (mostly IS200, 6–9 kb away). Its benzene dioxygenase / xylene monooxygenase / meta-cleavage operon on contig_0000378 (metabat_9342) is the most induced aromatic cluster, l2fc +6.5 in 2024-25 and +4.5 in 2025-26.
- **The exception is plasmid-borne alkane hydroxylases in Sphingorhabdus.** It carries 5.3% of the 2024-25 induction, 72.5% of it mobile:
  - **contig_0009679 (63.5 kb):** geNomad plasmid 0.98 with 16 hallmarks and a complete F-type conjugation system (TraA–TraG, MobF relaxase, TraD). It carries an IS3/transposase-flanked cassette: alcohol dehydrogenase, ferredoxin, CYP153 ("linalool 8-monooxygenase"), a TetR regulator, a transposase, alkB, an AraC (VirS) regulator and another transposase.
    - RNA in 2024-25: 1–3 CPM before oiling, 67 CPM at DAY07 (two days after diesel), peak 190 CPM at DAY12 (about 60× control).
    - RNA in 2025-26: 10–36 CPM in the oiled pool against about 1 in the control.
    - DNA: depth in the oiled pool rises from 0–0.5× before oiling to 3–5× by DAY07–10 (2024-25).
  - **contig_0014557 (47 kb):** a second F-type conjugative plasmid. It carries CYP153 + ferredoxin + acyl-CoA dehydrogenases (a β-oxidation module) beside an integrase/IS200/transposase cluster, and is induced about 20× in 2024-25.
- **Other catabolic plasmids:**
  - Mf105b01 (Parvibaculaceae, lorbin_068): a 132 kb plasmid with a P-type/virB T4SS, CYP153 and an intradiol dioxygenase. l2fc about +4.4 in 2024-25; copy ratio 1.0.
  - Yoonia: a 37 kb virB plasmid with ortho-cleavage genes, induced in 2024-25.
  - Pseudosulfitobacter: a 67 kb plasmid-like contig (score 0.97, no hallmark) with alkB, rubredoxin, RHO and both meta and ortho genes. It is the 2nd-largest plasmid contribution in both years (oiled 56–104 CPM vs control 23–33).
  - Pseudosulfitobacter: a separate 35 kb virB plasmid with alkB and alkG.
  - Rhodococcus: a 90 kb plasmid with CYP153.
- **Alcanivorax (2025-26):** its CYP153 on contig_0013262 (plasmid score 0.81, loose call) is 100% "mobile" but only 1.5% of the induction.

**Answer.** Community-level diesel degradation in both years is driven by chromosomally encoded pathways, overwhelmingly Cycloclasticus. HGT shows up in two ways:
- **Standing mobile pool:** IS-flanked alkB/CYP153 cassettes on F-type conjugative plasmids in alkane-degrading Sphingorhabdus (and Rhodococcus/Mf105b01), plus alkB and aromatic blocks shared across phyla. These respond fast, within 2 days in 2024-25, but are a minor share of total transcription.
- **Older acquisitions** inside Cycloclasticus. Per the cyclo fork, its dmp/tmo island proteins are 86–98% aa identical to C. zancles against a genome AAI of 82%. These are now fixed in the chromosome.

## Caveats
- **Pseudoreplication.** One pool per treatment. The 2024-25 EXP PRE01 sample already looks "oiled" (e.g. contig_0014557 at 26 CPM).
- **Annotation-based categories.** "RHO dioxygenase" includes non-hydrocarbon Rieske oxygenases. Bakta calls SAR11 Rieske proteins "biphenyl 2,3-dioxygenase". These make Pelagibacter 14.7% of the 2024-25 "induction", but the same genes are equally high before oiling and higher in the 2025-26 control, so the signal is pool composition.
- **Taxonomy.** Host taxonomy is Kaiju for 75% of HC contigs. Short contigs carrying a recently transferred gene can be mislabelled. The depth-correlation and flank tests are there to guard against this.
- **Coverage.** DNA depth comes from the CMO nanopore samples only to DAY10/11 (2024-25) and DAY6 (2025-26). With 10M RNA pairs per sample, most HC genes get 0 reads.

## Files (hgt/)
| file | contents |
|---|---|
| `hgt.json` | website data: counts, by_category, mobile_vs_chrom, mobile_share_by_category, mobile_genes, catabolic_mges (+ RNA/DNA courses), mobile_gene_courses, hgt_network (nodes/edges), hgt_flanks, sharing_funnel, induction_by_genus, plasmid_maps |
| `hc_genes.tsv`, `hc_context.tsv`, `hc_host.tsv`, `hc_expr_summary.tsv` | gene set, mobile flags, host taxonomy, expression |
| `shared_pairs.tsv`, `hgt_candidates.tsv`, `hgt_candidates_flanks.tsv` | sharing funnel |
| `catabolic_mges.tsv`, `mobile_hc_genes.tsv`, `mobile_vs_chrom_summary.tsv`, `mobile_share_by_category.tsv`, `induction_by_genus.tsv` | summary tables |
| `scripts/` | s1_select → s2_expr → s3_context → s4_host → s5_shared → s5b_depthcorr → s5c_candidates → s5d_flanks → s6_mobile_expr → s7_mge_view → s8_json |

The Grex side (extract.sh: proteins, regions, mmseqs) lives in ~/scratch/project_rnaseq/explore/hgt.
