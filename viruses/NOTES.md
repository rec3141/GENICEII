# Viruses in the CMO diesel mesocosms (2026-11-06)

Question: did the stress responses (phage-shock psp, envelope stress; Metals & toxicity page) come with a viral bloom after dosing?

## Viral set
- Inputs: geNomad virus_summary (272,410 contigs) and CheckV quality_summary, from mag_marine/mge on Grex.
- Confident set: virus_score ≥ 0.95, n_hallmarks ≥ 2, CheckV viral genes > host genes. That gives 44,809 sequences (91 geNomad provirus entries extracted as coordinate regions).
  - vOTUs: 42,192, clustered at ≥ 95% ANI over ≥ 85% of the shorter sequence (skani triangle -i, greedy longest-first; scripts/votu.py).
  - CheckV: 255 high-quality and 11 complete sequences; most are fragments.
  - Classes: Caudoviricetes 43,431; Megaviricetes 952 (giant viruses, kept separate); Maveriviricetes 138; Malgrandaviricetes 108; Faserviricetes 14.
- Some Megaviricetes contigs are binned into lorbin_356_sub and lorbin_507. lorbin_507, the "unclassified 5 Mb" bin whose islands GI_417/GI_3110 looked like prophages, holds 41 confident viral contigs, so it is partly a giant-virus/viral bin.

## Abundance (DNA)
- Nanopore, all 59 CMO pool libraries (mapping/depths.txt, corrected labels via common.nanopore_samples):
  - 2024-25: viral share of bases 5–14%.
  - 2025-26: 22–24% in both pools before dosing, falling to 13–16% after dosing in both. No diesel-specific rise.
- Illumina (first 2 M pairs of all QC-passing pool libraries mapped to the vOTU representatives): scripts/vmap.sh → map/counts → out/vdna_illumina.tsv.

## Activity (RNA)
- Viral CDS on confident non-provirus contigs, CPM of all CDS read pairs (scripts/vrna.py, vdiv.py).
- 2025-26 seawater medians:
  - diesel pool 46,950 before dosing → 10,970 after;
  - control 51,990 → 26,620;
  - active vOTUs (≥ 1 CPM): diesel 3,179 → 873, control 2,888 → 2,458;
  - Shannon: diesel 6.45 → 5.22, control 6.12 → 5.92.
- 2024-25: the diesel pool is higher than the control after dosing (11.6k vs 5.8k), but its single pre-dosing sample was 22.5k.
- Phages of bloomers (host = most co-varying MAG by Spearman across 130 libraries; scripts/vhost.py):
  - contig_0040828 (25.6 kb) ~ Cycloclasticus lorbin_722_sub, ρ 0.9. It goes from about 8–51 CPM before dosing to 935–1,229 at days 17–20, then 83 at day 82 while the host is still about 36k TPM. contig_0069512 (Schitoviridae) behaves the same way.
  - contig_0010334 (62 kb) ~ GCA-2708415 comebin_1134 (and Casp-alpha2), ρ 0.9. 2024-25 diesel pool only: 1,255 CPM before dosing, about 4,900 at day 8, about 900 at day 20. The host peaks at day 14. About ten more vOTUs follow the same host and are diesel-pool only.
- Giant viruses: in 2025-26 Megaviricetes go from about 5,800 to 250 CPM by day 17 in the diesel pool, against 4,700 at day 19 in the control. This is consistent with the algal decline.
- Virophage contig_0049668 (Priklausovirales, 22.8 kb, high quality): about 6,900 CPM in the 2025-26 control pool on days 3–5, absent from the diesel pool.
- Phage-shock: UBA9145 vamb_4149 has no geNomad virus or provirus and no Zot/pIV/filamentous-phage genes. Faserviricetes (Inoviridae) transcription is about 0 in the diesel pool. psp is not explained by phage.

## Prophages (scripts/vprophage.py)
- 441 geNomad proviruses, 164 on MAG contigs. Gene coordinates come from Grex annotation.tsv (reconciled tags), translated to local tags through locus_tag_map.
- The top diesel-pool "inductions" are not lytic:
  - Sulfurimonas lorbin_167 contig_0018983: host EPS/sugar genes at the provirus edge.
  - Unbinned contig_0001119: ParB/RepB/TniQ, plasmid- or Tn7-like.

## Virus-concentrate libraries (method test)
- VIRUS_NEG (no additive), PEG, MGCL2, PH and their combinations, 2025-26.
- Probable source: control-pool seawater from day 166 (5 May 2026, "Fernando's protocol", acidification) per the CMO 2025 log. This is inferred, not recorded on the libraries.
- Results:
  - viral transcripts 34–43k CPM, against 17–30k in control-pool seawater from days 49–81;
  - median Bray–Curtis 0.15 among treatments vs 0.42 to seawater;
  - untreated is as good as any additive;
  - enriched: virophage contig_0049668 (3,122 vs 159 CPM) and Phycodnaviridae/Megaviricetes.
- Kept out of every pool comparison.

## Synthesis
- No community-wide viral bloom: viral DNA, transcription and diversity fall after dosing, most in the diesel pool.
- Specific phages bloom with their oil-responding hosts: Cycloclasticus in 2025-26, the GCA-2708415/Casp-alpha2 group in 2024-25. This is consistent with host-density-driven ("kill-the-winner") dynamics. That interpretation is a hypothesis.
- The psp/envelope-stress signature is better explained by hydrocarbons than by phage.

## Illumina DNA (added after mapping; 90 QC-passing pool libraries, 2 M pairs each; out/vdna_illumina.tsv)
- Viral share of reads in late seawater:
  - 2024-25, days 25–70: diesel 8.7% vs control 16.7%;
  - 2025-26, days 12–15 (the last Illumina day): diesel 11.7% vs control 25.3%.
- vOTU richness, 2025-26 seawater, days 12–15: diesel about 8,100 vs control about 14,500.
- Cycloclasticus phage contig_0040828 DNA, diesel pool only: 0–8 reads per 4 M before dosing → 574 (day 11) and 468 (day 13).
- 2024-25 phage contig_0010334 DNA: already 1,100–3,350 reads in the diesel pool before dosing, with no rise. The RNA peak was activity.
- Virophage contig_0049668 DNA, control pool only: 7,418 (day 3) → 16,044 (day 12) and 15,906 (day 14) reads per 4 M; 0 in the diesel pool. This is the only clear viral bloom, and it is in the control pool.

## Bin domains (2026-11-07)
vhost and vprophage take their host labels from `special_label`, using the Tiara-based `data/bin_domain.tsv` (see the README).
