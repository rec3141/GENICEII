# Response guilds (2026-11-07)

Site page `#guilds` (Genomes & ice). Code:
- `scripts/guilds.py`: the analysis (permutations cached in `out/fingerprint_rows.pkl`; delete the cache to recompute).
- `scripts/guilds_site.py`: builds `site/guilds.json`.
- Input: the full-depth within-genome PC1–3 scores from `states/` (150 prokaryote MAGs, 108 core libraries).

## Methods
- **Layer 1, driver fingerprints**
  - Multivariate (RDA) R² of each genome's PC1–3 on five driver blocks:
    - winter;
    - ice (sea ice vs seawater, plus a top/bottom score);
    - diesel (diesel-pool indicator plus an "after dosing" indicator);
    - time (day, z-scored per winter);
    - light (log modelled under-ice PAR).
  - Commonality analysis gives each driver's unique R². "Seasonal" is the part shared among winter, time and light; "other shared" is the remainder.
  - Permutation tests (199 per genome) shuffle each block within strata that hold the other factors fixed:
    - winter: within medium;
    - ice: within winter;
    - diesel and light: within winter × medium;
    - time: within winter × medium × pool.
  - BH correction across genomes, per driver.
  - Embedding: UMAP (umap-learn 0.5) of the six fractions.
  - Clusters: sklearn HDBSCAN called every genome noise, so a Gaussian mixture is used (k = 5, silhouette 0.11; soft gradient). Clusters are named by the driver most over-represented relative to all genomes.
  - Single-winter genomes: the winter fraction is set to the median for the map only, and flagged.
- **Layer 2, co-response guilds**
  - MOFA+ is not installed (no installs), so the factor model is in numpy: PCA with missing values by EM (iterative SVD imputation) on the stacked genome scores. The matrix is 108 samples × (3 PCs × 150 genomes), and each genome block is scaled to unit total variance.
  - K = the largest number of factors before split-half stability drops below 0.7. Stability is the mean over 20 fixed random genome splits, comparing factor sample scores by best absolute correlation.

    | K | 1 | 2 | 3 | 4 | 5 |
    |---|---|---|---|---|---|
    | Stability | 0.95 | 0.88 | 0.80 | 0.71 | 0.66 |
    | Variance explained | | | | 72% | |

    That gives K = 4.
  - Genome loading on a factor = R² of its PC1–3 on that factor's sample scores. PC signs are arbitrary per genome, so magnitudes are used.
  - Guild = the factor with the highest z-scored loading. Clustering the loadings directly only split the two largest factors.

## Findings (full-depth RNA)
- **Drivers explain a median 67% of a genome's PC1–3 variance.**

  | Driver | Median unique share | Genomes with q < 0.05 |
  |---|---|---|
  | Ice vs water | 24% | 134 of 150 |
  | Diesel pool | 10% | 94 of 138 |
  | Winter | 8% | 65 of 87 |
  | Time | 2.3% | 38 of 150 |
  | Light | 2.1% | 32 of 150 |
  | Shared seasonal | 2.2% | – |

- **Fingerprint clusters:**

  | Cluster | Genomes |
  |---|---|
  | Diesel-pool responders | 57 |
  | Winter-switchers | 51 |
  | Time drifters | 16 |
  | Seasonal (shared drivers) | 15 |
  | Light followers, weakly structured | 11 |

- **Factors** (Spearman r with drivers; signs are arbitrary):

  | Factor | Tracks | Genomes in guild |
  |---|---|---|
  | F1 | winter (+0.87), ice (−0.43) | 36 |
  | F2 | ice vs water (−0.75), after dosing (−0.50), day (−0.44) | 35 |
  | F3 | day (+0.45), under-ice light (−0.39): season progression and darkening | 35 |
  | F4 | weak mix (diesel pool +0.31, ice −0.31, day −0.30) | 44 |

- **Key genomes:**
  - Cycloclasticus lorbin_184: diesel-pool responder (unique diesel R² 0.38), guild F4.
  - Cycloclasticus metabat_9342: winter-switcher, F1.
  - 50-400-T64 comebin_2576: ice-dominated, F2.
  - Sulfurimonas semibin_088: F3.
  - Planktomarina lorbin_118: diesel pool 0.25.
  - comebin_4318_sub: winter 0.56.
- **Same guild, different phylum** (cosine ≥ 0.97), candidates for coupling:
  - UBA9145 vamb_4149 with Verrucomicrobiota vamb_2202_sub (F2);
  - Actomonas metabat_1081 (Verrucomicrobiota) with 50-400-T64 comebin_2576 (F2);
  - Colwellia semibin_042 with Bacteroidota lorbin_317 (F1);
  - Actinomycetota SHLQ01 lorbin_111 with Porticoccus semibin_004 (F1).

## Caveats
- One pool per treatment, so "diesel pool" is a pool effect. In 2024-25 the pools were filled two weeks apart.
- Winter, time and light are collinear. The commonality split separates the unique and joint parts, but a small unique light fraction does not mean light is unimportant.
- 63 genomes are single-winter.
- Strain switches inside a MAG look like state changes.
- Clusters are soft: read the maps as gradients.
