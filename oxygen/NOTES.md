# Oxygen / low-O2 expression index (2026-11-07)

Site section `#oxygen` (site/oxygen.js, site/oxygen.json). There is no dissolved-O2 data: `do_pct` is empty.

## Code
- **scripts/oxygen.py**
  - Builds 20 O2-responsive gene groups from the corrected KOs:
    - terminal oxidases: cbb3, bd (high affinity, HA); aa3, bo3 (low affinity, LA)
    - paired enzymes: hemN/hemF, nrdDG/nrdAB
    - anaerobic markers: nar, nap, hox, group-4 hydrogenases, pfl, ldh, adhE
    - regulators fnr, arcAB, and globins
  - Community indices are log2 ratios of CPM.
  - Within-genome index:
    - Uses the 94 genomes with ≥300 HA+LA reads (≥30 of each type).
    - Per sample: log2(HA/LA) reads, needing ≥20 reads.
    - Each genome is centred on its own mean, then averaged weighted by sqrt(reads).
  - Composition vs regulation decomposition of the binned HA fraction:
    - fixed per-genome fractions × that sample's genome mix = composition;
    - the remainder = regulation.
  - Contrasts: paired ice − water and top − bottom, diesel − control, and BACI.
  - Covariates: Spearman correlations with ice physics (ice/phys/cores.tsv, nearest core ±2 d), water temperature, roof and mixers, and Cycloclasticus TPM.
  - Outputs:
    - out/group_cpm.tsv
    - out/sample_index.tsv
    - out/genome_shift.tsv
    - out/genome_ha.tsv
    - out/contrasts.tsv
    - out/covariates.tsv
    - out/within.pkl
- **scripts/oxygen_site.py** writes site/oxygen.json, which adds HA/LA carriers in 2025-26, Cycloclasticus lorbin_184 per-genome shares and the oil-snow comparison.

## Findings
**2025-26 seawater: the "low-O2" signal is composition.**
- Community HA/LA rises: diesel − control +1.42 log2 (p < 0.001), BACI +1.53.
- The genome mix alone predicts +1.15 of the difference (BACI +1.25).
- Regulation is −0.01 (BACI −0.09). The within-genome index is −0.12 (p = 0.50; BACI −0.13).
- Cycloclasticus lorbin_184 supplies 365 of 806 CPM of HA oxidase transcripts in the diesel pool. Its own cbb3/aa3 ratio is 0.40 (diesel) vs 0.40 (control) vs 0.31 (before dosing).
- Its hox share is 2.1‰ after dosing vs 3.35‰ before; narGHI is <1‰.
- So the cbb3/bd rise on the sulfur page, and the anaerobic-marker rise (+1.6 BACI), follow the bloom; they don't show a physiological switch.

**2024-25**
- Within-genome index: diesel − control +0.44 (p = 0.004).
- But BACI is −0.13: the single EXP pre-dosing sample was already higher, and the pools were filled 2 weeks apart. Most likely pool history.

**Season**
- In 2025-26 the within-genome index drifts up over the season in both pools (ρ 0.53 control, 0.68 diesel water). Seasonal, not diesel.

**Ice vs water**
- Community HA/LA: +0.50 in ice (22/24 pairs, p < 0.001). Mostly composition (+0.38, 21/24); Colwellia-type ice genera are cbb3/bd-rich.
- Within-genome: +0.10 (14/24, p = 0.37).
- Top vs bottom ice, within-genome excluding Cycloclasticus: +0.12 (9/13, p = 0.04), weak.
- Anaerobic markers are lower in ice (−0.58, p = 0.004).

**Covariates**
- The within-genome index correlates negatively with ice thickness in ice samples (ρ −0.59, n = 19) and with top-ice temperature (water, ρ −0.39), and positively with top brine salinity (ρ 0.48).
- All of these move together with season, so they can't be separated from time.
- With the roof open, the index is lower by 0.22 (p = 0.017; confounded with date).

**Paired enzymes**
- No consistent diesel signal: nrdDG/nrdAB is −0.97 diesel − control but +2.39 BACI, because the EXP pool started low. hemN/hemF is −0.10.

**Oil snow (CMO_EXP_DAY04_MOS)**
- More fermentation genes than same-day water and ice: pflB+pflA 84 vs 19/25 CPM, adhE 51 vs 0/23, nrdDG 11 vs 1/10.
- Oxidase indices are unremarkable. Single library.

## Limits
- The index detects regulation of the oxidase mix. It doesn't measure O2.
- It misses anoxic micro-niches (droplets, aggregates, brine) that hold few cells.
- About 78% of oxidase transcripts are binned.

## Recommendation
Install an under-ice O2 optode logger in each pool, and measure O2 in melted core sections and oil snow.

## Bin domains (2026-11-07)
Labels come from `special_label`, using the Tiara-based `data/bin_domain.tsv` (see the README). Eukaryotic and organelle bins are relabelled, and they are excluded wherever a bacterial genome is assumed. The three eukaryotic bins are comebin_4527, comebin_3937 and lorbin_494_sub.
