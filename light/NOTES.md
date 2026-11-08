# Light (2026-11-07)

Site section `#light`, first item under the Biogeochemistry menu (site/light.js, site/light.json). Code:
- **scripts/par_model.py**: modelled daily light at CMO Churchill (58.76° N, 94.07° W; local time UTC−6), written to out/par_daily.tsv for 2024-11-01 to 2025-03-15 and 2025-11-01 to 2026-03-15.
  - Solar position from the NOAA/Spencer equations, integrated in 1-minute steps.
  - Outputs: day length, top-of-atmosphere insolation and clear-sky GHI (Haurwitz: 1098 cosZ exp(−0.057/cosZ) W m⁻²).
  - PAR = 0.45 × GHI × 4.57 µmol J⁻¹, in mol photons m⁻² d⁻¹.
  - This is a clear-sky value above the ice and roof, so an upper bound: no cloud, snow, ice, roof or oil attenuation. No PAR was measured (the metadata `par` column is empty).
- **scripts/light.py**: 17 gene groups from KOs plus Bakta product regexes.
  - Oxygenic genes are restricted to phototroph hosts: plastid MAGs, Cyanobacteriota MAGs, or Kaiju-Cyanobacteria unbinned contigs, which are mostly algal plastids because Kaiju has no eukaryotes.
  - AAP genes are restricted to non-phototroph hosts, since bchG/bchP/acsF have chlorophyll-synthesis homologs.
  - Rhodopsins are split into proton pumps and heliorhodopsins.
  - Outputs:
    - per-sample CPM, summaries (diesel/control, BACI, ice, top/bottom, Spearman r with PAR within pool × season);
    - within-genome shares for AAP, rhodopsin and plastid genomes (out/genome_share.tsv, using genome CDS CPM from nitrogen/out/bin_cds_cpm.pkl);
    - per-sample date, PAR, day length and roof_open (out/sample_light.tsv).
- **scripts/light_site.py**: writes site/light.json, including the page text with numbers filled in from the data.
  - It also adds plastid photosynthesis share, metabat_1122 vs its states PC1 and N acquisition, and AAP vs rhodopsin partitioning by class.
  - It drops Kaiju-Synechococcus and plastid bins from the AAP genome list.
- **Dates**: samples are mapped by collection_date. Day 0 is 2024-11-28 and 2025-11-20 (the dominant collection_date − day offset); a few pre-dosing samples carry other offsets.

## Findings (seawater medians; descriptive, one pool per treatment)

**Modelled light.** Churchill has no polar night.
- At the solstice: 5.9 h of daylight, noon sun about 8° up, clear-sky PAR 2.5 mol m⁻² d⁻¹.
- By mid-February clear-sky PAR is about 13, roughly a 5-fold rise.
- Both experiments start near the darkest point, so PAR is confounded with experiment day.

**Algae.**
- Plastid photosystem II is −3.6 log₂ in the 2025-26 diesel pool vs control (BACI −3.4), following the drop in plastid TPM:

  | Plastid TPM | Diesel | Control |
  |---|---|---|
  | After dosing | 7,828 | 36,099 |
  | Before dosing | 46,373 | 39,963 |

- The photosynthesis share of plastid transcripts is unchanged: 813‰ vs 776‰.
- Within-genome BACI across 16 plastid/cyanobacterial genomes is 0.0.
- So the algae were lost, not shaded or switched off.
- Photosynthesis genes are lower in top ice than bottom ice (2 of 13 pairs higher at the top), consistent with bottom-ice algae.

**AAP.**
- Reaction-centre genes are +2.2 log₂ in 2025-26 (BACI +2.0), and LH1 +1.9.
- About 77% of reaction-centre transcripts are alphaproteobacterial (Sulfitobacter, RS62, Octadecabacter, Planktomarina).
- Within genomes the AAP share does not rise: median within-genome BACI −0.29, 2 of 8 genomes up. The community rise is composition.

**Rhodopsins.**
- Proton pumps carry about 10× more transcripts than AAP reaction centres.
- Hosts: alphaproteobacteria 36% (SAR11-type), actinobacteria 30%, Bacteroidota 17%.
- No diesel response: within-genome BACI +0.20, 11 of 20 genomes up.

**PAR correlations.** Within-genome correlations with modelled PAR are near zero (median r −0.13 to +0.13; plastids +0.31 in the diesel pool). There is no seasonal switch-on.

**Ice.**
- Higher in top ice than bottom ice:

  | Gene group | Pairs higher in top ice |
  |---|---|
  | BLUF photoreceptors | 13/13 |
  | PpsR/AppA | 13/13 |
  | photolyases | 12/13 |

- BLUF is also higher in ice than water in 19 of 24 pairs.
- Most AAP and rhodopsin genomes express less of their light genes in ice than in water (median within-genome log₂ ice − water: AAP −0.65, rhodopsins −0.28).

**metabat_1122 (Sulfitobacter).**
- AAP genes are 4.9‰ of its transcripts in water vs 2.8‰ in ice (p = 0.0009).
- In seawater:
  - AAP share follows the N-hungry end of its states PC1 (r −0.74);
  - it tracks its own N-acquisition share (r +0.46);
  - in 2025-26 it rises with PAR (r +0.46).
- Season, light and N hunger rise together, so these cannot be separated.

## Caveats
- Modelled PAR is a clear-sky upper bound and is confounded with experiment day.
- roof_open (54 libraries logged) is a pool-specific light modifier that the model does not include.
- Kaiju labels unbinned plastid genes as cyanobacteria.
- The `blh` and `phr` KOs are broad.

## Under-ice light model (2026-11-07)
- **Script:** `scripts/underice.py`. Per pool, clear-sky surface PAR is multiplied by a transmitted fraction: (1 − albedo) × exp(−15 m⁻¹ × snow − 1.5 m⁻¹ × ice).
  - Albedo: 0.85 with snow, 0.50 bare ice (> 2 cm), 0.07 open water.
  - No roof (as instructed), cloud or oil.
- **Inputs:** snow and ice from `ice/phys/cores.tsv` (CMO "Ice physical" sheets), linearly interpolated between core days and held flat outside them.
- **Outputs:** `underice_daily.tsv`, `underice_samples.tsv` and `underice_corr.tsv`.
- **Transmission:**

| Season | Pool | Transmitted |
|---|---|---|
| 2024-25 | both | 2–8% |
| 2025-26 | diesel | 90% → 0.6% |
| 2025-26 | control | 43% → 0.6% |

- **Correlations:** within-genome light-gene shares vs under-ice PAR have median r near 0.
- **metabat_1122 AAP share:** r −0.64 with under-ice PAR in the diesel pool (2025-26 water), −0.36 in the control; r with day +0.69. This reverses the "+0.46 with surface PAR" reading: AAP rises as the water darkens, which fits light repression of bacteriochlorophyll synthesis. It is still confounded with season and nitrogen hunger.
- **Roof state:** logged only on 2025-26 sampling days; drawn as pool-coloured vertical lines (solid = open, dotted = closed).

## Correlations with modelled under-ice light (rerun 2026-11-07)
- **Method:** `underice.py` now also writes `cat_under_corr.tsv` and `underice_partial.tsv`.
  - Spearman correlation with under-ice PAR and with transmittance, per pool.
  - A time-controlled partial: rank correlation after removing a linear trend in rank(day) within each pool, both pools pooled.
- **Result: nothing passes.**
  - Gene groups: smallest q is 0.22.
  - 180 within-genome shares: smallest q is 0.17.
  - Median partial r in 2025-26: AAP −0.14, rhodopsin +0.06, plastid photosynthesis −0.16.
- **Raw correlations are time trends.** Plastid photosystems vs under-ice PAR give r +0.74 to +0.81 in the 2025-26 diesel pool, because the algae collapsed as the ice thickened.
- **metabat_1122:** AAP vs under-ice PAR is r −0.64 raw but +0.04 with time held constant. The site text no longer claims a light response.
- **Framing (revised after user feedback):** "season" is not a cause, because microbes feel it through light, temperature, ice, brine and nutrients.
  - The time-controlled partial is a strict test for light effects beyond the shared seasonal trajectory. A null result does not rule light out.
  - metabat_1122 AAP tracks under-ice darkening (r −0.64) and not surface daylight (r +0.05). So its seasonal cue is darkening or another seasonal factor, not returning daylight.

## Bin domains (2026-11-07)
Labels come from `special_label`, using the Tiara-based `data/bin_domain.tsv` (see the README). Eukaryotic and organelle bins are relabelled, and they are excluded wherever a bacterial genome is assumed. The three eukaryotic bins are comebin_4527, comebin_3937 and lorbin_494_sub.
- **Full depth (2026-11-07):** partial correlations: smallest q is 0.07 for gene groups and 0.09 for genomes, still none < 0.05. metabat_1122 AAP vs under-ice r −0.68 (diesel), partial r +0.01.

## Day-18 exclusion (2026-11-08)
- In 2025-26 water, metabat_1122 AAP share is 5.9‰ in water vs 2.7‰ in ice (p < 0.0001; was 4.9 vs 2.9‰).
- Its correlation with N-acquisition share is r +0.33 (was +0.45). The time-controlled correlations still include no q < 0.05.
