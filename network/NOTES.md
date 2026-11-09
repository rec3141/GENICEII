# Genome network (2026-11-08)

Site page `#network` (Food web). Code in `scripts/`, run in order by `scripts/run_all.sh`. Small outputs in `out/`; large intermediate tables in `cache/` (untracked).

## Nodes
All node sets are read from files at run time; no MAG ids are hard-coded (except the six known pairs used only for validation, in `consensus.py`).

| Domain | Rule | Source |
|---|---|---|
| Prokaryote | DAS Tool MAG, `is_prokaryote` and not `is_plastid`, CheckM2 completeness ≥ 50% and contamination < 10% | `data/checkm2_dastool.tsv` (Grex `mag_marine/binning/checkm2/quality_report.tsv`, `dastool-*` rows) |
| Eukaryote | every EE## in `data/hull_mags.tsv`, any status; doubtful ones flagged | main checkout `data/` |
| Organelle | P##/M## groups (`data/organelle_groups.tsv`) and DAS Tool bins with `is_plastid`; groups over 1 Mbp are flagged multi-genome (a plastid genome is 0.1–0.4 Mbp; P01 is 5 Mbp) | |
| Virus | the Viruses page host table (`viruses/out/votu_hosts.tsv`) plus every Megaviricetes (giant virus) and Maveriviricetes (virophage) vOTU with mean viral-CDS CPM ≥ 5 over the 128 RNA libraries | `viruses/out/votu_rna.tsv` |

The site's own quality table (`data/summary.tsv`) only has DAS Tool SCG completeness. CheckM2 was copied from Grex (one `cat`, 2026-11-08).

**Contig assignment.** Each contig counts for one node only, in the order virus > organelle group > EE MAG > DAS Tool bin. For example, giant virus contig_0000610 is binned in plastid bin lorbin_567 and counts for the virus only.
- Five plastid bins have every contig inside P01: metabat_1343, lorbin_613, lorbin_663, lorbin_603 and metabat_7332_sub. They are dropped as nodes. So the known pair metabat_7332_sub → EE05 is tested as P01 → EE05.

Candidates: 528 prokaryotes (224 HQ, 304 MQ), 26 EE (after PR #18), 1 organelle group plus 27 plastid bins (22 with contigs of their own), and 63 vOTUs (26 giant viruses, 1 virophage, 36 phages). EE flags come from `status` and PR #18's `genome_status`. Doubtful: EE01, EE08 and EE16. Multi-genome: EE05, EE30 and P01.

## Abundance
- **RNA (since 2026-11-08: outlier-resistant profiles).** Every node's RNA profile is the shared robust profile of `euks/mags/scripts/robustexpr.py` (the Eukaryotes page metric):
  - per contig, reads per kb per million contig-assigned read pairs, 128 Illumina RNA libraries;
  - contigs detected in ≥ 25% of libraries, each scaled to mean 1;
  - per-library 10% trimmed mean × the node's median contig level.
  - The methods need counts, so the profile becomes read-pair pseudo-counts: profile × node kb × library size / 10⁶. Nodes with < 5 such contigs (every vOTU, small bins) fall back to the length-weighted mean, which returns their summed read pairs exactly.
  - It is applied to every node, prokaryotes included, so the layer uses one metric (cheap: one pass over the count matrix). The old summed counts are in `cache/rna_counts_summed.tsv`; `out/nodes_all.tsv` has `rna_contigs_used`, `rna_method`, `rna_thin` (< 20 contigs) and `rna_rho_robust_summed`.
- **DNA.** Nanopore coverage (depth × length) on the node's contigs from `euks/mags/dna/depths.txt.gz` (jgi depth table of the mag_marine mapping). Library size = coverage summed over all contigs.
  - Methods get pseudo-counts in kb of aligned bases.
  - **Why 44 DNA samples, not 72.** The depth table holds 276 libraries (the whole co-assembly mapping), and 72 of them are CMO runs listed in `dna/nanopore_samples.tsv`. Of those 72, 8 are source-site runs, 3 enrichments and 2 blanks, and 59 are pool runs (CTL/EXP). One pool run is a marine-oil-snow sample (CMO_EXP_DAY04_MOS, neither seawater nor ice), leaving 58 runs. Resequenced runs of the same sample are summed, giving 44 samples. Four pool runs in `nanopore_samples.tsv` (both PRE01s, a CTL_PRE02 and the EXP day-6 filtrate) are not in the depth table.
- **Prevalence rule, per layer** (`netcommon.py`):
  - detected in ≥ 25% of the layer's libraries and mean ≥ 1 per million;
  - EE MAGs, giant viruses, the virophage and organelle groups are always kept if detected in ≥ 10%.

  Result:

  | Layer | Libraries | Nodes | Prokaryote | Eukaryote | Virus | Organelle |
  |---|---|---|---|---|---|---|
  | DNA | 44 | 494 | 396 | 26 | 59 | 13 |
  | RNA | 128 | 590 | 483 | 26 | 63 | 18 |
- **Meta-variables** (FlashWeave): winter, diesel pool, control pool, seawater, sea ice, top ice, bottom ice (0/1), experiment day, and modelled under-ice PAR (`light/out/underice_*.tsv`; the DNA samples take the daily value of their winter, pool and day). Missing day/PAR (non-pool RNA libraries) are set to the winter median.
  - Water temperature is logged for 40 of 128 RNA libraries and none of the ice samples, so it is left out.

## Methods
1. **FlashWeave 0.19.2**: heterogeneous + sensitive mode (FlashWeaveHE-S), max_k 3, alpha 0.01, meta-variables as nodes, 6 workers. Defaults otherwise, including the convergence heuristic.
   - n_obs_min was set automatically to 20 non-zero samples: 91 DNA nodes and 1 RNA node fall below it and get no edges.
   - FlashWeave recommends HE mode for thousands of samples, and we have 44 and 128. It was used as requested. Read FlashWeave's DNA edges as weak.
2. **SPRING**: Python port in `scripts/spring.py` (see below), StARS with 20 subsamples of 80% and threshold 0.1, default lambda path. Edge sign = symmetrised MB coefficient.
3. **SPIEC-EASI** (R SpiecEasi, `scripts/run_spieceasi.R`): mb, 20 lambdas down to 0.01 × lambda.max, StARS with 20 subsamples (threshold 0.05), symBeta "ave".
4. **Residual Spearman** (`scripts/run_spearman.py`): centred log-ratio of (per million + 1); OLS on winter × pool × medium cells plus a linear day trend per winter; Spearman of the residuals; BH over all pairs; q < 0.05.
   - It is dense (see the table below). Spline time terms (cubic per winter or per winter × pool) only reduce RNA edges from 74k to 58k, so the residual structure is real covariation, not a missing time term.

**Consensus.** An edge needs:
- support from ≥ 2 methods with the same sign;
- no method with the opposite sign;
- at least one supporting method that accounts for the design (FlashWeave with meta-variables, or residual Spearman).

SPRING and SPIEC-EASI, like most compositional network methods, do not condition on winter, pool or medium, so their agreement alone can be a shared response to those drivers.

## SPRING port (`scripts/spring.py`)
Sources:
- SPRING (Yoon, Gaynanova & Müller 2019, Front. Genet. 10:516; github.com/GraceYoon/SPRING `R/SPRING.R`, `R/helpers.R`);
- latentcor (github.com/mingzehuang/latentcor `R/latentcor.R`, `R/internal.R`; the TT bridge also in mixedCCA `R/bridge.R`, Yoon, Carroll & Gaynanova 2020);
- `Matrix::nearPD`, `huge::huge.mb` and `pulsar::pulsar` / `stars.stability` / `.starsind` (read from the installed R packages).

Steps:
- **mclr** as in R.
- **Kendall tau-a** from sign products.
- **TT bridge:** two 4-variate normal CDFs, scipy's Genz.
- **Inversion:** a cached grid over r × z1 × z2 (39 × 39 zero fractions 0–0.95, 145 r values; `cache/bridge_tt_grid.npz`, about 5 min to build), bilinear in z and linear inverse in r.
  - Pairs with |tau| > 0.9 × bound go to root finding (latentcor's "approx" switch), as do pairs with a zero fraction above 0.95.
- **nearPD:** Higham with Dykstra, then the do2eigen step and nu = 0.001.
- **MB lasso:** coordinate descent on the correlation matrix, active set, warm starts, sym "or".
- **StARS:** cummax and `.starsind`.
  - The path is computed lambda by lambda over all subsamples and stops after the instability first reaches the threshold. Later lambdas cannot change pulsar's choice, so the result is the same.

Deviations: numpy RNG for subsamples, so the subsamples differ from R's; the grid replaces latentcor's stored interpolant.

Validation (`scripts/validate_spring.py` → `out/spring_validation.tsv`; R reference `scripts/spring_ref.R` with latentcor installed in `/data/scratch/grid-envs/Rlib` and SPRING sourced from a checkout).
- **Against R on the same simulated data** (n 120, p 30):
  - mclr and Kendall tau-a agree to 1e-14;
  - latent correlation ('original') agrees to 3e-5;
  - the grid matches R 'original' to 2e-4 and R 'approx' to 3e-3 (latentcor's own interpolant is the less exact of the two);
  - on the same correlation matrix the lambda path differs in 28 edge calls over 20 lambdas, all near-zero coefficients at small lambda (huge's looser tolerance);
  - StARS picks the same lambda (4/20); the refitted graph has 19 edges in Python and 18 in R, sharing 18 (Jaccard 0.95);
  - StARS instability differs by up to 0.012 (different subsamples).
- **Bridge grid against the exact bridge:** median r error 9e-5, max 1e-3. At z = 0 the bridge reproduces (2/π) asin r to 3e-5.
- **Simulation** (truncated Gaussian copula, compositional):

  | n | p | Precision | Recall | Sign agreement | Top-\|Pearson clr\| precision at the same edge count |
  |---|---|---|---|---|---|
  | 100 | 40 | 0.68 | 0.61 | 100% | 0.45 |
  | 200 | 60 | 0.66 | 0.61 | 100% | 0.54 |

- **Run time:** about 15 minutes per layer on a loaded 16-core machine, mostly the latent correlations of the 20 subsamples.

## Layout
t-SNE per layer (`scripts/layout.py`):
- nodes are points and libraries are features;
- values are the fourth root of the node's share of the library, with no other scaling;
- sklearn, PCA init, seed 20261108, perplexity = round(√nodes) (DNA 22, RNA 24).

## Re-running after EE changes
EE ids, statuses and taxonomy are read at run time from `data/hull_mags.tsv`, `data/contig_hull_mag.tsv.gz` and `data/ee_taxonomy.tsv`, so `run_all.sh` follows merges and relabels. This branch is rebased on `euk-robust-expr` (PR #20: robust RNA profiles, EE31 merged into EE03) and was rebuilt on its tables. Only the known-pair list in `consensus.py` names EE ids.

## Results (2026-11-08, robust RNA profiles; 25 EE MAGs after EE31 merged into EE03)
| Layer | Nodes (prok / euk / virus / organelle) | FlashWeave | SPRING | SPIEC-EASI | Residual Spearman (q < 0.05) | Consensus | Sign conflicts |
|---|---|---|---|---|---|---|---|
| DNA (44 samples) | 493 (396 / 25 / 59 / 13) | 740 | 5,694 (λ 0.053) | 3,611 (λ 0.143) | 3,834 | 1,253 | 8 |
| RNA (128 libraries) | 588 (482 / 25 / 63 / 18) | 1,855 | 10,051 (λ 0.042) | 6,424 (λ 0.088) | 69,711 | 8,067 | 214 |

- **Key-order bug fixed (affects the earlier results).** R (SPIEC-EASI) wrote each pair with the ids sorted case-insensitively, Python by code point, so a pair such as `contig_0000362`–`EE04` was `contig…, EE04` in SPIEC-EASI and `EE04, contig…` in the other three. Such pairs (every EE##/P## with a lower-case id) were counted as two different edges, so SPIEC-EASI never agreed with the other methods on them. `consensus.py` and `network_site.py` now put every method's pairs through `edge_key`. Consensus grew from 7,992 to 8,067 (RNA) and 1,240 to 1,253 (DNA). The earlier PR #19 tables had the same bug.
- **Robust RNA profiles:** Spearman of each node's robust vs summed RNA: prokaryotes median 0.976 (min 0.46), vOTUs 1.0 (single contigs fall back to the sum), EE MAGs median 0.93 (EE01 −0.13, EE03 0.68). Thin profiles (< 20 contigs): EE08, EE10, EE12, EE36 (flagged "thin RNA profile" on the page).
- **Residual Spearman** is still dense (69,711 RNA pairs at q < 0.05, 40% of pairs). Its edges keep ρ and q; the page uses them through the permissivity slider. Spearman-only edges are exported at |ρ| ≥ 0.4 (43,758 RNA edges below that are left out of `site/network_rna.json`; all are in `out/edges_spearman_rna.tsv`).

## Full networks and the permissivity slider (`network_site.py`, `site/network.js`)
- **Export:** `site/network_{dna,rna}.json` hold the nodes and the full edge sets of all four methods as column arrays: sign, masks of agreeing / opposing methods, FlashWeave |weight|, SPRING and SPIEC-EASI stability and |β|, Spearman |ρ| and −log₁₀ q. Score columns are sparse (one value per edge the method reports). `site/network.json` keeps the methods, the steps, summaries, known pairs and the partner tables.
  - Sizes: network.json 89 kB, network_dna.json 344 kB, network_rna.json 843 kB.
- **Steps** (`STEPS`; nested, each step is the union of its own clause and all stricter ones). "Strong" means SPRING/SPIEC-EASI stability ≥ 0.9, Spearman |ρ| ≥ 0.6 and q < 0.001; FlashWeave always counts. A step needs k of the switched-on methods, none with the opposite sign; design clauses also need FlashWeave or Spearman among the supporters. The method checkboxes still work (k is capped at the number switched on), as do the grazer, virus–host and plastid–host views.

| Step | Rule | DNA links | RNA links |
|---|---|---|---|
| 1 | all four, strong | 15 | 283 |
| 2 | three, strong | 151 | 870 |
| 3 | three | 564 | 3,114 |
| 4 | three, or two strong (one design-aware) | 646 | 3,344 |
| 5 | two, one design-aware (consensus, default) | 1,253 | 8,067 |
| 6 | any two | 2,296 | 8,759 |
| 7 | any two, or one strong | 3,203 | 11,548 |
| 8 | any one | 10,901 | 31,082 |

- **Known pairs recovered** (`out/permissivity_steps.tsv`; not used for tuning):
  - RNA: P01–EE05 from step 1; contig_0000362–EE04 from step 3; contig_0000564–EE07 from step 5; contig_0000610–EE22 and contig_0000623–EE07 from step 7; virophage contig_0049668–EE07 at no step. Its only support is a Spearman edge below the |ρ| 0.4 export floor.
  - DNA: P01–EE05 from step 2, contig_0000623–EE07 from step 3, the others at no step.
  - So recall in RNA rises from 1/6 at steps 1–2 to 2/6 at step 3, 3/6 at the consensus and 5/6 at steps 7–8. Step 7 costs 43% more links than the consensus.
