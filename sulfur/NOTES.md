# Sulfur & hydrogen (2026-11-06)

Site section `#sulfur` (site/sulfur.js, site/sulfur.json). Code:
- scripts/sulfur.py: 45 gene groups from corrected KOs, plus a Bakta-product regex for dsrEFH/dsrD.
  - Per-sample CPM of CDS reads goes to out/cat_cpm.tsv; per-gene CPM to out/gene_cpm.pkl; genes to out/genes.tsv.
  - Hosts go to out/cat_hosts.tsv; summaries (diesel vs control, 2025-26 BACI, ice vs water, top vs bottom) to out/cat_summary.tsv; dsr/apr/qmo hosts to out/dsr_genes.tsv.
- scripts/sulfur_site.py: adds Kaiju genera for unbinned genes, the dsr direction table and the Sulfurimonas block, and writes site/sulfur.json.
  - It drops groups with no genes, such as dsrD (no KO, no product hit).
- Gene IDs are the old BAKTA_EXTRA tags, consistent with the local counts.

## Findings (seawater medians after dosing; log2 diesel/control; descriptive, one pool per treatment)
- **Sulfurimonas hydrogen oxidation and CO2 fixation:**
  - The quinone-reactive NiFe hydrogenase (hydAB, Campylobacterota type) is +4.5 (2024-25) and +3.1 (2025-26), with 2025-26 BACI +3.1.
  - Diesel-pool transcripts split 48% Sulfurimonas semibin_088 and 50% unbinned Sulfurovum-like contigs (Kaiju).
  - aclAB is 96% Sulfurimonas (+4.0 / +4.0, BACI +2.1). porABCD is 99% Sulfurimonas (+4.8 / +3.3, BACI +2.0).
  - fccAB sulfide dehydrogenase is +2.7 in 2025-26 (BACI +1.8): 58% semibin_088, 15% Cycloclasticus.
- **Pool history vs diesel:**
  - 2025-26: both pools had similar hydAB before dosing (about 2,950 vs 2,870 CPM). The diesel pool kept it while the control fell, so this winter shows a real diesel effect (BACI q=0.023 for semibin_088 at MAG level).
  - 2024-25: the single diesel pre-dosing sample already held 903 CPM vs 15–69 in the control (after the day-5 relabel). Here the contrast is mostly pool history; MAG BACI is −0.3 to −1.3.
- **Cycloclasticus is a hydrogen/sulfur actor too:**
  - It carries 93% of the bidirectional Hox hydrogenase (+5.7 in 2025-26, BACI +3.8).
  - Also: membrane nitrate reductase narGH (+3.2, BACI +3.9; Cycloclasticus sp019090295 and sp001566825), ETHE1 (+2.5), mddA (+3.0/+3.1; 63% Cycloclasticus), soxCD and dsrEFH.
- **Sulfate reduction is a trace activity:**
  - Reducer phyla carry ≤7% of dsrAB, ≤1% of aprAB and ≤6% of qmo transcripts.
  - dsrAB comes from Gammaproteobacteria sulfur oxidisers (Thiohalobacter- and Thermochromatium-like contigs, GCA-001735895, UBA9214); aprAB mostly from SAR11.
  - The only SRB MAG is Desulforhopalus comebin_3581_sub.
- **Electron acceptors:**
  - High-affinity oxidases rise in 2025-26 diesel: cbb3 +2.2 (BACI +2.0), bd +1.1.
  - Denitrification steps do not: norBC −2.8, nosZ +0.1, napAB −0.5.
  - This reads as an O2-limitation signature with Cycloclasticus nitrate respiration, not Sulfurimonas denitrification.
  - **Revised 2026-11-07 (oxygen/):** within-genome oxidase ratios do not shift; the cbb3/bd and narGH rise is the Cycloclasticus bloom (composition), so this is not evidence of O2 limitation.
  - The logs have no usable dissolved-oxygen data: one 2024-25 pre-dosing value, 105%.
- **Sulfurimonas genomes:**
  - semibin_088 (1.7 Mbp, 92%): sqr, fccAB, hyaAB/hydAB, aclAB, PFOR, cbb3 and aa3, but no sox and no nap/nar. It is a microaerophilic H2/sulfide oxidiser.
  - lorbin_167 (3.1 Mbp, 84%): sox, nap, nir, nor, nos and group-4 hydrogenase. The denitrifying generalist.
  - comebin_1140_sub resembles semibin_088.
  - All three are ice-enriched in RNA (17/24 pairs) and DNA (semibin_088 18/19).
  - Pool FST 0.35 (semibin_088) and 0.23 (comebin_1140_sub), against ice/water FST 0.03–0.04.
- **DMSP/methanethiol:** the dmdBCD pathway (mostly Sulfitobacter) is +1.0 / BACI +1.9 in 2025-26. DMSP lyases are mixed.
- **RuBisCO:** −1.9 / −3.0, from algal plastids, consistent with the algae declining in the diesel pool.
- **Hypotheses (labelled as such on the page):**
  - H2/sulfide sources: anoxic micro-niches in oil, oil snow or brine; fermenter and Hox H2; organosulfur in the diesel and from DMSP.
  - Low-O2 brine fits Sulfurimonas's ice preference.
