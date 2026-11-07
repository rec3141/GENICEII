# Physiological states (2026-11-07)

Site section `#states` (site/states.js, site/states.json), in the "Diesel response" nav group.

## Code
- **scripts/states.py** (about 15 min; loads the full gene matrix):
  - per-MAG within-genome profiles and PCA;
  - discrete states (GMM BIC, k-means silhouette);
  - PERMANOVA on the top 5 PCs, permuting within medium;
  - BACI-style shift ratio;
  - ice vs water;
  - changepoints;
  - function correlations;
  - growth index;
  - binomial thinning check.
- **scripts/states_post.py**: what the top-3-PC space aligns with (season, pool, medium R² and ARI), plus PC and growth trends with day per pool.
- **scripts/states_site.py**: builds site/states.json, including the data-driven "What the states show" bullets.
- **Outputs** in out/: mag_states.tsv, functions.tsv, changepoints.tsv, rarefy.tsv, assoc.tsv, trends.tsv; scores, loadings, growth and depth as pkl files.

## Method
- **Sample and MAG selection:**
  - Core samples only: CTL/EXP, water or ice.
  - A MAG is used in samples where it has at least 3,000 CDS reads, and only if that holds in at least 20 samples. 85 of 911 MAGs pass.
- **Genes:** kept if detected in at least 90% of the used samples with a mean of at least 5 reads; top 3,000 by mean.
- **Transform:** log2 within-genome CPM with a depth-scaled prior (0.5 reads at mean depth), centred per sample, then PCA.
  - The first run used a fixed pseudocount and genes detected in at least 50% of samples. PC1 then tracked depth strongly (|r| up to 0.9).
- **Depth checks:**
  - With the current settings, |r(PC1, log depth)| > 0.6 in 6 of 85 MAGs.
  - Binomial thinning to each MAG's minimum depth gives PC1 agreement |r| < 0.7 in 3 of 85 MAGs.
  - For the key MAGs, agreement is 0.72–0.98.
- **GMM BIC is not used in the text.** It favours more than one state in 79 of 85 MAGs, which is unreliable with n = 20–100 in 3 dimensions.
  - Instead, the page reports which design factor the PC space aligns with.

## Findings
- **Ice vs water** is the main axis in 58 of 85 MAGs (winter 14, pool 13).
  - Ice/water PERMANOVA after dosing is significant in 53 of 60 MAGs (2024-25) and 45 of 47 (2025-26).
- **Diesel vs control after dosing, 2025-26:** only 20 MAGs are testable, because blooming MAGs have too few reads in the control. 17 of 20 have q < 0.05; median R² is 0.19.
  - Top: Planktomarina lorbin_118 (R² 0.62, shift 2.1×; PC1 = ribosomal +0.86 vs cold shock −0.82) and RS62 metabat_1104 (0.58, 1.7×).
- **Sulfurimonas semibin_088:**
  - R² is 0.45, but the shift ratio is 0.36, so this is a pool (strain) difference, not diesel.
  - PC1 is ribosomal +0.92 vs motility −0.94.
  - The growth index falls through 2025-26 in both pools: 98→38‰ in diesel, 115→27‰ in control.
- **Cycloclasticus lorbin_184:**
  - PC1 rises with day in the diesel pool (r +0.66).
  - Along PC1, hydrocarbon genes (−0.79) and storage (−0.58) fall while chaperones rise (+0.57).
  - The growth index falls 67→44‰.
- **Cycloclasticus metabat_9342:** PC1 is a growth axis (ribosomal +0.90, cold shock +0.77) against a foraging axis (motility −0.77, transporters −0.64, storage −0.65, N acquisition −0.59).
- **50-400-T64 comebin_2576:** the growth index rises 69→81‰ over the 2025-26 diesel season. Its main axis pairs chaperones, ROS and cold-shock genes against transporters and N acquisition.
- **Background heterotrophs slow in diesel only, 2025-26:**

  | MAG | Diesel growth index | Control growth index | Other stats |
  |---|---|---|---|
  | Sulfitobacter metabat_5376 | 137→96‰ | 160→157‰ | shift 2.3×, R² 0.12, q 0.001 |
  | UBA9145 vamb_4149 | 81→52‰ | 74→80‰ | psp vs PC1 r −0.77 |

- **Not testable:** Nitrosopumilus metabat_9833 does not pass the read threshold.

## Caveats
- One pool per treatment, so diesel cannot be separated from other pool differences.
- 2024-25 pools were filled two weeks apart.
- Different strains in different pools or winters look like different states.
- The growth index (ribosomal-protein share) is a rough proxy.

## Depth-aware refit (2026-11-07)
- **Problem:** large PCA outliers were mostly shallow libraries. 62 of 74 outliers (robust z > 6) had fewer than 10k reads on the genome; for comebin_2576, all three were its shallowest libraries (~3.6–4.2k reads vs a median of 69k).
- **Fix:** PCA axes are now fitted only on samples with ≥ 10,000 reads on the genome (FIT_READS; all samples if fewer than 15 qualify). Samples with 3,000–10,000 reads are projected onto those axes. The site draws them small and faded.
- **Effect:**
  - Outliers 74 → 35 (shallow 62 → 23); comebin_2576 has none.
  - Ice vs water is now the main axis for 65 of 85 MAGs (was 58).
  - 2025-26 pool differences: 19 of 20 (was 17).
  - metabat_8566: ice R² 0.59 → 0.65, depth r 0.52 → 0.38.
  - lorbin_184: no steady PC1 trend with day (r −0.08), so the text now rests on the growth index.
- **Remaining deep outliers are coherent groups:** late-season 2025-26 diesel ice (metabat_8566, vamb_856_sub) and 2024-25 day-55 diesel samples (Colwellia lorbin_095, vamb_500_sub).
- **Tooltip:** now spells out the ice layer (top/bottom section, etc.).
