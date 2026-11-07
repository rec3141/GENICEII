# Nitrogen (2026-11-07)

Site section `#nitrogen` (site/nitrogen.js, site/nitrogen.json). Code:
- **scripts/nitrogen.py:** 30 gene groups from the corrected KOs.
  - narGH in Nitrospira/Nitrospina hosts counts as nxrAB.
  - amoABC counts only in nitrifier hosts; other CuMMO copies are kept separate.
  - Outputs:
    - per-sample CPM: out/cat_cpm.tsv
    - summaries: out/cat_summary.tsv
    - per-genome CDS CPM: out/bin_cds_cpm.pkl (cached)
    - host-normalised N-acquisition share: out/host_nstress.tsv
- **scripts/nitrogen_site.py:** writes site/nitrogen.json, with:
  - hosts (MAG or Kaiju genus)
  - time courses
  - nstress
  - the nitrifier/CuMMO genomes block
- **No nutrient data.** The CMO metadata has no NH4 or NO3 values. Filtrate was collected from 2025-26 onward but not analysed.

## Findings (seawater medians after dosing; descriptive, one pool per treatment)

**N-scavenging programme in the 2025-26 diesel pool** (CPM diesel vs control; BACI in log2)

| Gene group | Diesel | Control |
|---|---|---|
| nitrate/nitrite transporters (BACI +4.5) | 2,197 | 70 |
| nasA/narB | 802 | 22 |
| nirBD | 865 | 27 |
| amtB | 3,633 | 627 |
| glnB/K | 1,216 | 159 |
| glnA | 4,284 | 1,066 |
| urtA–E | 1,580 | 114 |
| cynS | 97 | 0 |

- Carried by Cycloclasticus, Sulfurimonas semibin_088 and Sulfitobacter.
- Host-normalised:
  - Cycloclasticus N-acquisition share doubles (log2 +1.0).
  - Sulfitobacter genomes rise by log2 +1.5 to +2.5.

**2024-25 is pool history.** The single EXP pre-dosing sample already had 3,816 vs 85 CPM of nitrate transporters, so the BACI is negative.

**Nitrate is reduced but not denitrified.**
- Cycloclasticus narGHI is +3.2 (BACI +4.5) and hcp +2.6.
- norBC is −2.8, nosZ stays low, and nirK/nirS are flat.

**Other N pathways**
- N2 fixation is below 10 CPM (Arcobacter-like) and shows no diesel response.
- Anammox genes are absent. nrfA is at trace levels.

**Nitrifiers (2024-25 only)**
- Nitrosopumilus metabat_9833 plus Nitrosomonadaceae comebin_369 total about 60 TPM in the diesel pool vs 190 in the control (pre-dosing: 190 vs 380).
- The per-genome amo share is unchanged (about 300‰).
- MAG BACI is −0.6 and −1.6 (not significant).
- Nitrifiers are nearly absent in 2025-26.

**CuMMO in hydrocarbon degraders**
- Two genomes carry pmo/amo-like operons:

| Genome | TPM | BACI | q |
|---|---|---|---|
| Cycloclasticus metabat_9342 | 8,140 | +3.7 | 0.009 |
| Porticoccaceae 50-400-T64 comebin_2576 | 35,915 (≈3.6% of transcripts) | +3.6 | 0.02 |

- These are likely short-chain alkane monooxygenases (Rubin-Blum et al. 2017, doi:10.1038/nmicrobiol.2017.93).
- Pitfall: GTDB puts Cycloclasticus in Methylococcales, so a naive order-based "methanotroph" split misassigns it.

**Top vs bottom ice**
- Nitrate transporters are higher in top ice in 0 of 13 pairs; urea transporters in 1 of 13.
- One hypothesis is that brine concentrates nutrients.
