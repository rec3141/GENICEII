# CMO RNA-seq exploration (2026-10-23)

130 BGI metatranscriptomes (10M pairs/sample) mapped to the marine Flye coassembly
(illumina_rna run in Grex ~/scratch/project_rnaseq/marine). Metadata from
/data/dev/genice-sample-metadata. EXP = oiled pool, diesel both years: 2024-25 added DAY05 16:30; 2025-26 ~136 L injected DAY01
(crude oil was only used for isolates/small incubations). CTL = unoiled pool. One pool each, so time points are
pseudoreplicates and p-values are descriptive.

## Website
Private artifact "Churchill Diesel Mesocosms": https://claude.ai/artifact/VctR8TksDcnKt4WySD6t49
- Source: `site/` (index.html, app.js with D3, data.json, cyclo.json, hgt.json).
- Rebuild data with `scripts/build_site.py`, which needs `pathways.py`, `pathways_site.py`, `ko_cpm.py`, `dna.py` and `hc.py` run first.
- Narrative text is in `scripts/site_text.py`.
- Publish from the scratchpad mirror `cmo-site/`.

## Layout
- `data/`: expression tables pulled from Grex, plus `genemap.tsv.gz` (mtx row → gene, contig, KO, DAS Tool bin),
  GTDB-Tk r226, DAS Tool summary, Kaiju contigs and the KEGG KO list.
- `scripts/`: `common.py` (loaders, derived design: season/pool/medium/phase), `load_genes.py`, `hc.py` (hydrocarbon
  genes), `who.py` (which contigs/MAGs carry them), `mags.py` (MAG oiled vs control), `baci.py` (pool × period),
  `ordination.py` (Aitchison PCA + PERMANOVA), `kos.py` (KOs consistent across seasons), `figs.py`.
- `out/`: result tables (`pathways/` = module completeness and expression per MAG). `fig/`: static PNGs from the first pass.
- `dna/`: per-MAG nanopore DNA depth (Grex job explore/dna/bindepth.sh) and the CMO run→sample map.
- `canthyd/`: CANT-HYD hydrocarbon HMM hits (trusted cutoffs) for all 8.7 M proteins.
- `cyclo/`: Cycloclasticus genomics deep dive (fork). See cyclo/NOTES.md.
- `hgt/`: mobile context and cross-taxon sharing of hydrocarbon genes (fork). See hgt/NOTES.md.
- `audit/`: independent check of the site's claims (AUDIT.md). The site text was corrected from it.
- `mgscreen/`: k-mer screen of the 173 CMO Illumina metagenomes (local /matika reads, first 400k R1 reads) against 15 key MAGs and 2 plasmids. Tool: `scripts/kscreen.py`; aggregation: `scripts/mg.py` → `out/mg_kmer.tsv`.
- `ISSUE_ko_misalignment.md`: root cause and fix options for the pipeline KO misalignment (mag_analysis basic vs extra Bakta).

## Main findings (updated 2026-10-23 evening, corrected KOs)
See the website (claude.ai artifact "Churchill Diesel Mesocosms") and `site/` for the full version.
1. Season dominates (PERMANOVA R² 0.50), then pool (R² 0.27–0.45 within season); ice vs water is minor.
2. Cycloclasticus dominates aromatic degradation in both winters: 83–99% of tmo/tod/nah/extradiol/catE transcripts in the diesel pool.
   - 2024-25: metabat_9342 (98.7% ANI to an Arctic sea-ice MAG), rising after day 20 to 19–51% of mapped transcripts.
   - 2025-26: lorbin_184 + sibling sub-bins, already resident in both pools before dosing (nanopore DNA 0.3–0.6% of binned in CTL, 1.0–1.8% in EXP), dominant from day 2.
   Source-water sites had no Cycloclasticus DNA.
3. In 2025-26 seawater the diesel/control expression ratios are: tmoA-F 43×, todC 37×, extradiol DO 26×, nahAc 16×, catE 16×, alkB-like 3–7×. BACI-corrected (vs own pre-dosing baseline) the aromatic responses are 11–14×.
   CANT-HYD (trusted cutoffs) gives canonical AlkB to Pseudosulfitobacter and CYP153 to unbinned Alcanivorax/Sphingorhabdus.
   The Cycloclasticus alkB-like genes fall below the AlkB cutoff.
4. A Sulfurimonas H2/sulfide chemolithoautotrophy program rises in the diesel pool: porA, NiFe hydrogenase, aclA (rTCA), fccB, all ≥95% Sulfurimonas.
   In 2025-26 it was already high pre-dosing in both pools. dsrAB stays low.
5. Cu membrane monooxygenase (pmoA-like) of Porticoccaceae 50-400-T64 reaches 30–40 CPM in the 2025-26 diesel pool (≈0 in control). Propane MO (prmA) of unbinned Pseudonocardia is up 6–10× in 2024-25, but only about 1.5× after correcting for the diesel pool's pre-dosing level.
   Animal signalling genes, HSP90 and PSII PsbZ fall in every comparison.
6. Pools differed before dosing. Mislabelled libraries were corrected on 2026-10-26 (see `swapcheck/` and the archive's sources/sample_corrections.tsv):
   - The 2024-25 day-5 extract tubes were swapped (RNA, Illumina and 7 nanopore runs relabelled).
   - The Illumina PRE02 pair (A03/A04) was swapped.
   - Illumina A05 and C08 are wrong-pool libraries without a partner, and 8 Illumina libraries failed; all of these are excluded.
   - After the rerun, 2024-25 pool R² is 0.454 (was 0.416). The diesel pool's 2024-25 RNA baseline is now one sample (EXP_PRE01).
7. KO tables from the pipeline (ko_tpm/cog/cazy, and the public RNA viz) are misaligned. The mag_marine eggNOG/kofam run used the Sep-24 basic Bakta IDs.
   Here KOs are remapped by coordinates (`data/ko_remap.tsv.gz`, `data/ko_cpm.tsv`).

8. HGT is mostly a side show. 2% of confident hydrocarbon genes are on MGEs, and mobile-associated genes carry about 9% of the diesel-induced transcription.
   The standout is a Sphingorhabdus 63.5 kb F-type conjugative plasmid (contig_0009679) with an IS3-flanked CYP153+alkB cassette, about 60× induced.
   26 flank-confirmed cross-genus transfers were found, 7 of them expressed.
9. Marine oil snow samples (2024-25 day 62/63) are hotspots: tmo 3,300–5,500 CPM, up to 3× day 55–60 seawater; Cycloclasticus 32–51% of transcripts.
10. Illumina metagenome k-mer screen (173 libraries, `out/mg_kmer.tsv`, values = % of sampled read k-mers unique to the target):
   - metabat_9342: ≈0.01% (noise floor) everywhere until day 25, then 10–14% of diesel-pool k-mers at days 55–63 (water and ice); control 0.01% all season. June 2025 diesel-pool surface: 0.03% (residual).
   - lorbin_184: diesel pool 0.13% pre-dosing → 0.87% (SW, days 12–25); control 0.05–0.08%.
   - Sphingorhabdus plasmid contig_0009679: 10–20× more in diesel-pool DNA at days 5–25 of 2024-25 (≈0.02% vs 0.001%).
   - Sulfurimonas semibin_088: 2024-25 diesel pool ≈0.01% → 1–3% (days 5–25); 2025-26 it is 3–12% in both pools.
   - Source/incubation/MOS/SERF metagenomes show no Cycloclasticus above the floor, so the bloom seed is not traceable at 400k reads. The crude-oil incubations favour Oleispira lorbin_369 (up to 0.13%) instead.

## QC notes
- In 2025-26 controls, 55–65% of assigned pairs land on rRNA loci that bbduk misses. Most are the euglenid 18S
  contig_1340804 (Phacus-like, about 13M pairs in CTL_SW). So "% rRNA" understates rRNA there, and CDS yield is about
  17% vs about 41% in EXP. Expression values are CDS-only, so TPM is unaffected; depth is lower.
- Unbinned contigs carry 36% (median) of mapped pairs.
- 17 Kaiju contig lineages in mag_marine hold peptide strings: taxon IDs 3242192/3242193 are missing from names.dmp.
- Placement by profile similarity: the 8 VIRUS_* libraries match 2025-26 CTL seawater, days 5–19 (r≈0.96), and are
  near-identical to each other. CMO_MYSTERY matches 2024-25 CTL sea-ice TOP, days 6–11. These are candidate dates for
  the metadata gaps.

## Ice vs water, top vs bottom (2026-10-24)
Code: `ice/scripts/` (phys.py, contrasts.py, rrna_contrasts.py, panel.py, divergence.py, ice_site.py → site/ice.json, rendered by site/ice.js).
Outputs: `ice/phys/` (cores.tsv, profiles.tsv from the "Ice physical" log sheets), `ice/out/` (paired contrasts), `ice/cold/` and `ice/popgen/` (fork workstreams).
Python with pandas on this machine: /data/miniforge3/bin/python3 (system python3 is 3.8 without pandas).
- Design: 24 ice–water pairs (same season/pool/day; 2024-25 TOP+BOT averaged, 2025-26 whole cores) and 13 top–bottom pairs (2024-25). The pairs are repeated measures in one pool per treatment, so treat sign counts as consistency only.
- Physics: top ice −16 to −2 °C, bottom −6.5 to −1.6 °C. Brine salinity at the top is 61–166 PSU vs about 29–107 at the bottom. Brine volume is 3–10% at the top vs 7–25% at the bottom, so the top is near or below the 5% permeability threshold. 2025-26 CTL day 12 temperatures are inverted in the log and were dropped.
- MAGs, ice vs water (CLR TPM): ice-enriched MAGs are Cognaticolwellia (4 species, all +0.4 to +1.7), Glaciecola, Pseudomonas_E, Nonlabens, Psychrobacter, Shewanella vesiculosa, Marinomonas primoryensis, Polaribacter vamb_3778, Maribacter arcticus and Flavobacterium. Water-enriched MAGs are Methylophilaceae BACL14, Pseudohongiellaceae UBA9145, Opitutaceae and Haliangiaceae.
  - Siblings split: Colwellia (5 bins) lean to water/bottom ice, while Cognaticolwellia leans to ice/top.
  - Sulfitobacter is mixed: most species are top-leaning, vamb_549 is bottom-leaning.
- Plastids: the "Gloeomargaritales" MAGs are 100–280 kb chloroplast genomes, identified by SILVA plastid 16S. Hosts include Attheya septentrionalis (an ice diatom, lorbin_567), Guillardia (lorbin_511), Costaria (kelp) and Thalassiosira. Bins are listed in data/plastid_bins.tsv, and common.plastid_label() relabels them on the site.
- KO panel, ice > water:
  - polyamine synthesis 24/24 pairs;
  - cspA 20/24 and DEAD-box helicases 20/24;
  - gas vesicles 21/24;
  - peroxiredoxins 17/24, ectoine 17/24 and osmolyte uptake 17/24;
  - methylotrophy is lower in ice.
- KO panel, top > bottom:
  - Dps/osmY 13/13;
  - branched-chain fatty-acid start (bkd) 12/13, matching membrane fluidity in colder top ice;
  - storage 13/13 and DMSP/mddA 13/13;
  - bottom ice has more RuBisCO (12/13) and light harvesting (11/13).
- rRNA taxa:
  - ice: Chlamydomonas (24/24), Cryothecomonas, Protaspis, Melosira and Nitzschia/Cylindrotheca (top), plus terrestrial deposits (Picea pollen 24/24, Trebouxia, Dothideales fungi);
  - water: Phaeocystis, Prymnesium, Katablepharis and Euglenozoa;
  - bottom ice: dinoflagellates, Ebria and cryptophytes.
- Ice diverges from its paired water by an Aitchison distance of about 17–19 (seawater replicates about 9), and the distance grows with ice thickness (2025-26 rho 0.73).

## Day-5 label correction rebuild (2026-10-26)
- Another session added cross-platform identity checks (swapcheck/, identity.tsv; common.corrections()) and the QC flags in the metadata. It reran hc, baci, mags, ordination, kos, pathways, dna, figs, mg and the ice scripts at 06:37–06:56, then republished.
- The cyclo/ and hgt/ pipelines read the old labels and were rerun at 13:07 (log: swapcheck/rerun_cyclo_hgt.log; pre-fix copies in swapcheck/pre_fix/cyclo_hgt/). cyclo/scripts/{dna,dna_rel}.py and hgt/scripts/s7_mge_view.py now use common.nanopore_samples(), which applies the nanopore corrections.
- Changes are small, because day 5 is pre-dosing.

## Genomic islands (2026-10-26)
Code: islands/scripts/islands.py (annotation + expression), islands_site.py (→ site/islands.json), site/islands.js (page).
Grex helper: ~/scratch/project_rnaseq/explore/islands/islands.sh (Bakta features on island contigs + island/contig GC), local copies in islands/grex/.
- IslandPath-DIMOB (mag_marine/mge/islandpath): 3,863 islands, 2,146 in 494 MAGs. 179 contain at least 2 core genes (ribosomal proteins, RNAP, EF-Tu, ATP synthase). These are codon-usage false positives, flagged `core_like`.
- Inside vs outside on the same contigs, enriched: integrase 9.6×, transposase 7.3×, R-M 4.2×, phage 4.2×, defence/TA 3.0×, metal resistance 1.3×.
- Depleted: oil genes 0.4× (51 of 10,775 confident genes, in 31 islands), transport 0.5×, EPS 0.5×.
- Oil islands:
  - Colwellia semibin_042 GI_348: a full dmpK–P phenol hydroxylase with MopR, next to an integrase and transposase. Silent.
  - Cycloclasticus metabat_5432 GI_2286: nahB next to IS66. 13× raw in the 2025-26 diesel pool, but its share of host transcripts is unchanged.
  - Parasphingorhabdus GI_1266: mhp meta-cleavage operon. 18× in the 2024-25 diesel pool, host share unchanged.
  - GI_2179: CYP153 + ferredoxin reductase + IS21, in an unclassified bin. Switched on within its host (0→3 CPM, 19× host per-gene average).
- Other islands responding to diesel within their host: Oleispira GI_734 (IS1182 + BON domain, 119 CPM, 39× host average). Cycloclasticus flagellin and type IV pilus islands are highly expressed (460–540 CPM). metabat_9342 has an ectoine island (GI_935).
- The lorbin_507 prophage islands (GI_417, GI_3110) are season-specific, not a diesel response: higher in diesel in 2024-25, much higher in the control in 2025-26.
- Other HGT markers in mag_marine:
  - geNomad plasmid/virus/provirus and IntegronFinder are fine (already used in hgt/).
  - MacSyFinder was run with --db-type unordered on the whole-assembly protein file, so its 131k "systems" have no locus information and cannot place conjugation systems.
  - DefenseFinder output is header-only (completed in 90 s on 13 M proteins); see ISSUE_mge_tools.md.

## Population genetics, ice vs water (2026-10-26)
Code: ice/popgen/scripts. The pipeline is map.sh (bwa mem vs all 911 MAGs, MAPQ≥20, ≥95% id, ≥80% aligned) → abundance.py → pileup_bin.py (bcftools mpileup) → popgen_bin.py → summarize.py.
- Samples: 55 libraries at 2 M pairs each. QC-flagged libraries are in samples_excluded_qc.tsv. Targets: 72 MAGs (targets.txt).
- Outputs: popgen_summary.tsv, paired_fst_vs_null.tsv, fst_outlier_genes.tsv, abundance_icewater_dna.tsv; site data in ice.json → popgen.
- No within-species ice/water structure: for 22 of 23 well-covered species, the median same-day ice-vs-water FST equals the consecutive-day water-vs-water FST (median excess −0.002). Top vs bottom ice also matches the baseline.
- The exception is MAG-120531 lorbin_317 (Flavobacteriaceae): FST 0.16 vs 0.08 baseline, and lower π in ice (0.0051 vs 0.0074).
- The pools hold different strains (ice, control vs diesel FST): Cognaticolwellia lorbin_135 0.37, Sulfurimonas semibin_088 0.35, Colwellia 0.27–0.28, UBA9145 0.17. Low-coverage MAGs reach 0.7 or more (Sulfitobacter vamb_549, UBA5078, BACL14 lorbin_054).
- Ecotypes sit at the species level. In DNA:
  - Cognaticolwellia lorbin_135 is higher in ice in 19/19 pairs; metabat_8244 18/19; Glaciecola 18/19. All three are top-enriched.
  - Colwellia vamb_500_sub and lorbin_095 lean to water and bottom ice.
  - Sulfitobacter species mostly lean to top ice; vamb_549 leans to the bottom.
- The 2025-26 pair FSTs (about 0.7) come from low coverage (<700 sites) and are not interpretable.

## Sulfur & hydrogen, Metals & toxicity pages (2026-11-06)
- Sulfur & hydrogen: sulfur/ (scripts/sulfur.py, sulfur_site.py; out/; NOTES.md) → site/sulfur.json + sulfur.js. Supersedes the "Beyond hydrocarbons" markers subsection on the pathways page, which is now a pointer. Overview finding 4 was rewritten: 2025-26 BACI is a real diesel effect; 2024-25 is mostly pool history.
- Metals & toxicity: metals/ (scripts/msets.py, mexpr.py, mmob.py, mbuild.py; grex/metal_ctx.sh; out/; NOTES.md) → site/metals.json + metals.js. No metal measurements exist; exposure is inferred from transcription.
- Bakta locus-tag fix check: audit/bakta_fix/NOTES.md. The site is unaffected (coordinate remap). Grex mag_marine extra files were reconciled on Nov 2 and the RNA tables regenerated on Nov 4. Local IDs are the old tags; translate through locus_tag_map.tsv.

## Nitrogen (2026-11-07)
`nitrogen/` (scripts, out, NOTES.md) feeds the site page `#nitrogen`. It covers:
- the diesel-driven N-scavenging response
- nitrate reduction without full denitrification
- the nitrifier decline
- copper monooxygenases in Cycloclasticus and Porticoccaceae 50-400-T64

See nitrogen/NOTES.md.

## 50-400-T64, oxygen, eukaryotes (2026-11-07)
- `porticocc/`: Porticoccaceae 50-400-T64 genome comebin_2576, the second-largest diesel responder of 2025-26.
  - Its CuMMO forms its own lineage, most likely an alkane monooxygenase.
  - Its toolkit is constitutive: diesel selected the organism rather than inducing its genes.
  - Site page `#porticocc`.
- `oxygen/`: low-O2 expression index, at community level and within each genome.
  - The diesel-pool shift towards high-affinity oxidases is composition (the Cycloclasticus bloom), not regulation.
  - This revises the O2-limitation reading on the sulfur page. Site page `#oxygen`.
- `euks/`: eukaryotes from 18S rRNA and plastid transcripts.
  - Algae collapse in the 2025-26 diesel pool and bacterivorous chrysophytes and choanoflagellates rise.
  - Giant viruses track their algal hosts; the control-pool virophage tracks Phaeocystis.
  - Site page `#euks`.
- Site nav is now grouped into dropdowns (NAV in site/app.js).

## Physiological states (2026-11-07)
`states/` holds per-MAG PCA of within-genome expression for 85 MAGs. Site page `#states`, under the Diesel response menu.
- Ice vs water is the dominant state axis.
- In 2025-26, 17 of 20 comparable MAGs differ between pools.
- Cycloclasticus lorbin_184 winds down as its bloom ages; 50-400-T64 keeps growing.
- Background Sulfitobacter and UBA9145 slow down in the diesel pool only.
- The Sulfurimonas pool difference is strain, not diesel.

See states/NOTES.md.

## Light (2026-11-07)
`light/` covers light use and feeds the site page `#light` (first item under Biogeochemistry). See light/NOTES.md.
- **Gene groups:** oxygenic photosynthesis, aerobic anoxygenic photosynthesis (AAP), rhodopsins, and light sensing/protection.
- **Modelled light:** clear-sky PAR for Churchill (`par_model.py`), since no PAR was measured.
- **Algae:** diesel removes the algae, but surviving plastids keep their photosynthesis share.
- **AAP:** its rise in the diesel pool is composition, not induction.
- **Rhodopsins:** dominant and flat.
- **Top ice:** expresses more light sensing.
- **Sulfitobacter metabat_1122:** AAP tracks its nitrogen-hungry water state.

## Bin domains (2026-11-07)
`scripts/bin_domains.py` → `data/bin_domain.tsv` classifies every bin from the contig-level Tiara calls in the microscape contig table. Rule, length-weighted: eukarya ≥ 50% → eukaryote; organelle (organelle/plastid/mitochondrion) ≥ 50% → organelle; otherwise prokaryote. Result: 885 prokaryote, 23 organelle, 3 eukaryote.

Helpers in `scripts/common.py`:

| Helper | Use |
|---|---|
| `bin_domain(b)` | the bin's call |
| `is_prokaryote(b)` | bacterial-genome filter |
| `is_plastid(b)` | `plastid_label` rule, or a Tiara organelle call |
| `special_label(b)` | plastid label, "Plastid/organelle (Tiara)" or "Eukaryote (Tiara)" |

Every page's label function now calls `special_label`. Analyses that assume a bacterial genome skip non-prokaryote bins:
- physiological states (per-MAG PCA and growth index)
- the oxygen within-genome index
- nitrogen host_nstress
- light AAP/rhodopsin hosts

The Eukaryotes page lists the eukaryotic nuclear bins. Rerun `bin_domains.py` whenever the bins change.
