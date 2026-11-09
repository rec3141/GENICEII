# Genome network (2026-11-08)

Site page `#network` (Food web). Code in `scripts/`, run in order by `scripts/run_all.sh`. Small outputs in `out/`; large intermediate tables in `cache/` (untracked).

## Nodes
All node sets are read from files at run time; no MAG ids are hard-coded (except the six known pairs used only for validation, in `consensus.py`).

| Domain | Rule | Source |
|---|---|---|
| Prokaryote | DAS Tool MAG, `is_prokaryote` and not `is_plastid`, CheckM2 completeness ≥ 30% and contamination < 20% (`MIN_COMPL`, `MAX_CONTAM`; since 2026-11-09, before ≥ 50% / < 10%), and expressed in the pools (below) | `data/checkm2_dastool.tsv` (Grex `mag_marine/binning/checkm2/quality_report.tsv`, `dastool-*` rows) |
| Eukaryote | every EE## in `data/hull_mags.tsv`, any status; doubtful ones flagged | main checkout `data/` |
| Organelle | P##/M## groups (`data/organelle_groups.tsv`) and DAS Tool bins with `is_plastid`; groups over 1 Mbp are flagged multi-genome (a plastid genome is 0.1–0.4 Mbp; P01 is 5 Mbp) | |
| Virus | the Viruses page host table (`viruses/out/votu_hosts.tsv`) plus every Megaviricetes (giant virus) and Maveriviricetes (virophage) vOTU with mean viral-CDS CPM ≥ 5 over the 128 RNA libraries | `viruses/out/votu_rna.tsv` |

The site's own quality table (`data/summary.tsv`) only has DAS Tool SCG completeness. CheckM2 was copied from Grex (one `cat`, 2026-11-08).

**Contig assignment.** Each contig counts for one node only, in the order virus > organelle group > EE MAG > DAS Tool bin. For example, giant virus contig_0000610 is binned in plastid bin lorbin_567 and counts for the virus only.
- Five plastid bins have every contig inside P01: metabat_1343, lorbin_613, lorbin_663, lorbin_603 and metabat_7332_sub. They are dropped as nodes. So the known pair metabat_7332_sub → EE05 is tested as P01 → EE05.

Candidates (2026-11-09): 703 prokaryotic bins pass 30/20 (94% of the bins' TPM, against 83% with 50/10); the pool expression rule keeps 496 of them (207 dropped, 149 of which passed the old 50/10 rule). 26 EE (after PR #18), 1 organelle group plus 27 plastid bins (22 with contigs of their own), and 63 vOTUs (26 giant viruses, 1 virophage, 36 phages). EE flags come from `status` and PR #18's `genome_status`. Doubtful: EE01, EE08 and EE16. Multi-genome: EE05, EE30 and P01.

**Prokaryote flags.** "fragment" when completeness < 50% (54 kept MAGs; drawn pale) and "contamination 10–20%" when contamination ≥ 10% (67; dotted ring); 5 are both. Shown in the hover and in `out/nodes_all.tsv` (`fragment`, `contam_10_20`). The high-contamination bins are mostly near-complete bins with a second strain, such as the Cycloclasticus bloom genomes lorbin_184 (94.7% / 10.2%) and metabat_5432 (100% / 12.8%).

**Pool expression rule** (`nodes.py`: `POOL_TPM_MEAN`, `POOL_TPM_MAX`, `POOL_DET`; both layers, before the per-layer prevalence rule). Pool libraries: `common.meta()` group CTL or EXP and medium SW or SI, both winters, before and after dosing (106 libraries); site, enrichment, MOS and other libraries do not count.
- Prokaryotic MAGs: bin TPM (`data/marine_coassembly.bin_tpm.tsv`) ≥ 10 on average over the pool libraries, OR ≥ 100 in at least one (keeps short bloomers).
- EE MAGs, organelle groups and plastid bins: robust RNA profile detected in ≥ 10% of the pool libraries (chosen over a percentile of the prokaryotes' values: simpler, and on the same profile the RNA layer uses). No EE MAG or plastid with contigs falls out.
- vOTUs: unchanged.
- Dropped nodes are in `out/nodes_dropped_pool_expression.tsv`; their contigs are not handed to another node.
- Kept, as expected: both Cycloclasticus bloom genomes, metabat_9342 and comebin_3787_sub (Cycloclasticus), comebin_4318_sub (MarineAlpha4-Bin1), comebin_1103_sub (Shewanella psychromarinicola), lorbin_206 (Pseudoalteromonas arctica), metabat_1243_sub (Hyphomonas), lorbin_516_sub (Yoonia). lorbin_783 (Moritella, 18% complete), metabat_8071_sub (Sulfitobacter, 18%) and lorbin_722_sub (Cycloclasticus, 23%) are below 30% completeness, so they are not candidates at all.

**GTDB labels fixed (2026-11-09).** `nodes.py` looked GTDB-Tk calls up by the plain binner name only; `*_sub` bins exist in `data/gtdbtk_taxonomy.tsv` only as `dastool-*` rows, so 182 of them were "unclassified" (e.g. comebin_1103_sub is Shewanella psychromarinicola). It now reads the `dastool-` row first. Where both rows exist they are identical.

**Possible same genome** (`network_site.py`: `SAME_RHO`, `SAME_N`). Two prokaryotic MAGs of the same GTDB genus (or, both without a GTDB genus, the same majority Kaiju genus) whose robust RNA profiles have Spearman ρ ≥ 0.95 over the ≥ 10 libraries where both are detected. Their links (either layer) are flagged; the page draws them dashed and can hide them. Pairs: `out/same_genome_pairs.tsv`.

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
5. **Presence Pearson** (`scripts/run_pearson.py`, added 2026-11-08):
   - for each pair, only the libraries where both nodes are > 0; Pearson r of log(share of the library), the same input as the other methods (`netcommon.load_layer`, `rel_abund`);
   - tested when co-present in ≥ 10 libraries (`PP_N_MIN`); t test with n − 2 df; BH over all tested pairs of the layer. Computed with matrix products (n = PᵀP, sums of log and log² masked by the partner's presence), checked against scipy `pearsonr` on sampled pairs.
   - **Edge call** (`netcommon.PP_R, PP_N, PP_Q`): **|r| ≥ 0.85, co-present in ≥ 20 libraries, q < 0.05.** Each edge keeps r, n and q (`out/edges_pearson_{layer}.tsv`). DNA 4,479 edges (3.7% of pairs), RNA 1,994 (1.2%).
   - Why not q alone: q is already BH-corrected, but with up to 128 libraries |r| ≈ 0.25 passes q < 0.05, so q < 0.05 alone gave 48,553 DNA (40% of pairs) and 104,840 RNA edges (61%). Rule chosen (2026-11-08, by the user's instruction "raise r or lower the p-value until it's sensible"): the least strict combination, raising |r| before lowering q, at which the count is at most SPRING's in both layers (DNA 5,694, RNA 10,051). |r| ≥ 0.8 passes RNA (4,019) but not DNA (7,470 at n ≥ 20, 6,692 even at q < 1e-5); |r| ≥ 0.85 with n ≥ 20 passes both. Lowering q hardly changes RNA (every |r| ≥ 0.7 pair already has q < 1e-5). An alternative with lower |r| would be |r| ≥ 0.75 with n ≥ 30 (DNA 4,445, RNA 6,962), which mainly removes rarer DNA genomes; not used.
   - Edges at every threshold (`out/pearson_thresholds.tsv`, all with q tested at n ≥ 10):

| n | \|r\| | DNA q<0.05 | DNA q<1e-3 | DNA q<1e-5 | RNA q<0.05 | RNA q<1e-3 | RNA q<1e-5 |
|---|---|---|---|---|---|---|---|
| ≥ 10 | ≥ 0.5 | 42,637 | 23,761 | 10,603 | 38,473 | 37,983 | 35,787 |
| ≥ 10 | ≥ 0.6 | 33,804 | 23,155 | 10,603 | 22,543 | 22,379 | 21,967 |
| ≥ 10 | ≥ 0.7 | 22,810 | 18,596 | 10,535 | 11,070 | 11,006 | 10,953 |
| ≥ 10 | ≥ 0.75 | 17,349 | 15,173 | 9,744 | 7,041 | 7,014 | 6,981 |
| ≥ 10 | ≥ 0.8 | 12,220 | 11,338 | 8,041 | 4,045 | 4,037 | 4,020 |
| ≥ 10 | ≥ 0.85 | 7,669 | 7,433 | 5,779 | 2,006 | 2,003 | 1,997 |
| ≥ 10 | ≥ 0.9 | 3,795 | 3,772 | 3,293 | 677 | 677 | 677 |
| ≥ 20 | ≥ 0.5 | 30,062 | 19,462 | 9,254 | 38,288 | 37,960 | 35,784 |
| ≥ 20 | ≥ 0.6 | 22,398 | 18,856 | 9,254 | 22,387 | 22,356 | 21,964 |
| ≥ 20 | ≥ 0.7 | 14,562 | 14,297 | 9,186 | 10,984 | 10,983 | 10,950 |
| ≥ 20 | ≥ 0.75 | 10,874 | 10,874 | 8,395 | 6,991 | 6,991 | 6,978 |
| ≥ 20 | ≥ 0.8 | 7,470 | 7,470 | 6,692 | 4,019 | 4,019 | 4,017 |
| ≥ 20 | ≥ 0.85 | 4,479 | 4,479 | 4,430 | 1,994 | 1,994 | 1,994 |
| ≥ 20 | ≥ 0.9 | 2,072 | 2,072 | 2,072 | 675 | 675 | 675 |
| ≥ 30 | ≥ 0.5 | 13,997 | 10,428 | 5,225 | 38,111 | 37,893 | 35,768 |
| ≥ 30 | ≥ 0.6 | 10,053 | 9,822 | 5,225 | 22,289 | 22,289 | 21,948 |
| ≥ 30 | ≥ 0.7 | 6,207 | 6,207 | 5,157 | 10,941 | 10,941 | 10,934 |
| ≥ 30 | ≥ 0.75 | 4,445 | 4,445 | 4,366 | 6,962 | 6,962 | 6,962 |
| ≥ 30 | ≥ 0.8 | 2,888 | 2,888 | 2,888 | 4,005 | 4,005 | 4,005 |
| ≥ 30 | ≥ 0.85 | 1,602 | 1,602 | 1,602 | 1,991 | 1,991 | 1,991 |
| ≥ 30 | ≥ 0.9 | 622 | 622 | 622 | 674 | 674 | 674 |
   - Not design-aware.
**Consensus.** An edge needs (five methods since 2026-11-08):
- support from ≥ 2 of the 5 methods with the same sign;
- no method with the opposite sign;
- at least one supporting method that accounts for the design (FlashWeave with meta-variables, or residual Spearman).

SPRING, SPIEC-EASI and the Presence Pearson do not condition on winter, pool or medium, so their agreement alone can be a shared response to those drivers.

**Effect of adding the Presence Pearson (same rule, strict edge call).** RNA consensus 8,067 → 8,665 (+598), DNA 1,253 → 1,490 (+237); no link of the four-method consensus is lost. Links supported only by residual Spearman + Presence Pearson: 584 of the RNA consensus (6.7%) and 97 of DNA. Sign conflicts: RNA 214 → 227, DNA 8 → 45. (With the first, q-only call, q < 0.05: RNA consensus 38,310, 79% of it Spearman + Presence Pearson only, 9,342 conflicts.)

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

## Layouts
- **t-SNE** per layer (`scripts/layout.py`): nodes are points and libraries are features; values are the fourth root of the node's share of the library, with no other scaling; sklearn, PCA init, seed 20261108, perplexity = round(√nodes) (DNA 22, RNA 24).
- **Force-directed**, one per method and one for the consensus, per layer (`scripts/forcelayout.py` → `forcelayout.R`):
  - edges: the method's own edges as the page exports them (after the floors), positive links only (a negative link is no reason to draw two genomes together);
  - weight (attraction) = the method's score, divided by its median: FlashWeave |weight|, SPRING and SPIEC-EASI StARS stability (floored at 0.02; SPRING keeps a few refit edges with stability 0), residual Spearman |ρ|, Presence Pearson |r|; consensus edges: number of supporting methods + mean percentile of their scores;
  - backbone: each genome keeps its 10 strongest positive links (a link stays when it is in the top 10 of either end). Without it the dense sets (residual Spearman: 28k positive RNA links) were one round ball;
  - igraph 2.1.4 `layout_with_fr` (R), started from the t-SNE map, 2000 iterations, no grid, seed 20261108: the same on every run;
  - genomes with no positive link in that layout go on rings outside it (radius 1.1 × the 99th percentile of the laid-out radius), ordered by domain then abundance (DNA FlashWeave: 116 genomes; DNA Spearman 33; Presence Pearson DNA 210, RNA 224; DNA consensus 22; RNA FlashWeave 5).
- Positions are stored per layer in `network_{layer}.json` (`layouts`), each scaled to −100…100, 1 decimal (RNA: 44 kB for 7 layouts). The page animates between layouts (750 ms; links are hidden during the move when more than 12,000 are drawn).
- What the force layouts show (Louvain on the layout backbones, for description only): consensus modularity 0.68 (DNA) and 0.57 (RNA). In RNA one module of 82 genomes holds 15 of the 25 EE MAGs and 11 of the 18 plastid nodes (with 50 prokaryotes and 6 vOTUs); in DNA one of 37 holds 12 EE MAGs. The other modules are almost all prokaryotes. Top positive consensus hubs: RNA comebin_4104_sub (Pararhizobium, 49 links), lorbin_111 (43), comebin_3924 (41); DNA lorbin_029 (Parasphingorhabdus, 19).

## Adding a layer (e.g. Illumina DNA)
`netcommon.LAYER_CFG` lists the layers: abundance table (nodes × libraries, in `cache/`), library table (`size` plus the design columns season, group, medium, day, layer, par_under), pseudo-count unit, page label and order. Every Python script loops over it, `run_spieceasi.R` and `run_flashweave.jl` run every `cache/*_X.tsv` that `export.py` writes, and the page builds its layer buttons and tables from `network.json`'s `layer_list`. A new layer needs: its two tables written by `nodes.py` (or by hand) and an entry in `LAYER_CFG`; prevalence is computed from the table when `nodes_all.tsv` lacks `{layer}_prev`. Built for Illumina DNA on 2026-11-09 (layer key `mgdna`, below).

## Re-running after EE changes
EE ids, statuses and taxonomy are read at run time from `data/hull_mags.tsv`, `data/contig_hull_mag.tsv.gz` and `data/ee_taxonomy.tsv`, so `run_all.sh` follows merges and relabels. This branch is rebased on `euk-robust-expr` (PR #20: robust RNA profiles, EE31 merged into EE03) and was rebuilt on its tables. Only the known-pair list in `consensus.py` names EE ids.

## Results (2026-11-08, robust RNA profiles; 25 EE MAGs after EE31 merged into EE03)
| Layer | Nodes (prok / euk / virus / organelle) | FlashWeave | SPRING | SPIEC-EASI | Residual Spearman (q < 0.05) | Presence Pearson (\|r\| ≥ 0.85, n ≥ 20, q < 0.05) | Consensus (4 methods) | Consensus (5 methods) | Sign conflicts (5) |
|---|---|---|---|---|---|---|---|---|---|
| DNA (44 samples) | 493 (396 / 25 / 59 / 13) | 740 | 5,694 (λ 0.053) | 3,611 (λ 0.143) | 3,834 | 4,479 | 1,253 | 1,490 | 45 |
| RNA (128 libraries) | 588 (482 / 25 / 63 / 18) | 1,855 | 10,051 (λ 0.042) | 6,424 (λ 0.088) | 69,711 | 1,994 | 8,067 | 8,665 | 227 |

Presence Pearson overlap (Jaccard; sign agreement), DNA / RNA: FlashWeave 0.079 / 0.150 (99.5% / 100%; 383 of FlashWeave's 740 DNA and 501 of its 1,855 RNA edges), SPRING 0.074 / 0.056 (99.6% / 99.8%), SPIEC-EASI 0.092 / 0.093 (99.7%), residual Spearman 0.040 / 0.021 (89% / 99%).

- **Key-order bug fixed (affects the earlier results).** R (SPIEC-EASI) wrote each pair with the ids sorted case-insensitively, Python by code point, so a pair such as `contig_0000362`–`EE04` was `contig…, EE04` in SPIEC-EASI and `EE04, contig…` in the other three. Such pairs (every EE##/P## with a lower-case id) were counted as two different edges, so SPIEC-EASI never agreed with the other methods on them. `consensus.py` and `network_site.py` now put every method's pairs through `edge_key`. Consensus grew from 7,992 to 8,067 (RNA) and 1,240 to 1,253 (DNA). The earlier PR #19 tables had the same bug.
- **Robust RNA profiles:** Spearman of each node's robust vs summed RNA: prokaryotes median 0.976 (min 0.46), vOTUs 1.0 (single contigs fall back to the sum), EE MAGs median 0.93 (EE01 −0.13, EE03 0.68). Thin profiles (< 20 contigs): EE08, EE10, EE12, EE36 (flagged "thin RNA profile" on the page).
- **Residual Spearman** is still dense (69,711 RNA pairs at q < 0.05, 40% of pairs). Its edges keep ρ and q; the page uses them through the permissivity slider. Spearman-only edges are exported at |ρ| ≥ 0.4 (43,758 RNA edges below that are left out of `site/network_rna.json`; all are in `out/edges_spearman_rna.tsv`).

## Full networks and the permissivity slider (`network_site.py`, `site/network.js`)
- **Export:** `site/network_{layer}.json` holds the nodes, the layouts and the edge sets of all five methods as column arrays: a, b; `f` packs the sign and the bit masks of agreeing, opposing and strong-scoring methods; sparse score columns (one value per edge the method reports): FlashWeave |weight|, SPRING and SPIEC-EASI stability, Spearman |ρ|, Presence Pearson |r| (×100). `site/network.json` keeps the methods (with their line-style score classes), layouts, layers, steps, summaries, known pairs and the partner tables.
  - Floor for single-method links: residual Spearman |ρ| ≥ 0.4 (43,488 RNA Spearman-only edges below it left out). The Presence Pearson needs no floor, its edge call is strict. All edges stay in `out/edges_*_{layer}.tsv`.
  - The page recomputes nothing about "strong"; it reads the bits, which `netcommon.strong_bits` sets on the exported (rounded) values.
  - Sizes: network.json 92 kB, network_dna.json 343 kB, network_rna.json 663 kB (was 843 kB with four methods; the edge flags are now packed into one column).
- **Strong** scores: SPRING/SPIEC-EASI stability ≥ 0.9; Spearman |ρ| ≥ 0.6 and q < 0.001; Presence Pearson |r| ≥ 0.9, q < 0.001 and co-present in ≥ 30 libraries (stricter than its edge call on all three; DNA 622, RNA 674 edges); FlashWeave always counts.
- **Steps** (`STEPS`, five methods; nested, each step is the union of its own clause and all stricter ones). A step needs k of the switched-on methods, none with the opposite sign. With five methods, three or four supporters can all be design-blind (SPRING, SPIEC-EASI, Presence Pearson), so every clause down to the consensus now also needs FlashWeave or the residual Spearman among its supporters; without that, the union would put design-blind links into step 5. `network_site.py` asserts that step 5 equals `consensus.py`'s consensus. The method checkboxes still work (k is capped at the number switched on). With some methods switched off, the sign is the one the switched-on methods give, so a method outvoted overall still shows its own links when it is the only one on (before, the 8 RNA FlashWeave links outvoted by other methods disappeared with FlashWeave alone).

| Step | Rule | DNA links | DNA known pairs | RNA links | RNA known pairs |
|---|---|---|---|---|---|
| 1 | four or five, strong | 54 | 0 | 377 | 1 |
| 2 | three strong, one design-aware | 185 | 1 | 939 | 1 |
| 3 | three, one design-aware | 717 | 2 | 3,454 | 2 |
| 4 | three, or two strong (one design-aware) | 808 | 2 | 3,670 | 2 |
| 5 | two, one design-aware (consensus, default) | 1,490 | 2 | 8,665 | 3 |
| 6 | any two | 2,986 | 2 | 9,400 | 3 |
| 7 | any two, or one strong | 4,046 | 2 | 12,163 | 5 |
| 8 | any one | 14,070 | 2 | 31,774 | 5 |

- **Known pairs** (`out/permissivity_steps.tsv`; not used for tuning; out of 6):
  - RNA: P01–EE05 from step 1; contig_0000362–EE04 from step 3; contig_0000564–EE07 from step 5; contig_0000610–EE22 and contig_0000623–EE07 from step 7; virophage contig_0049668–EE07 at no step. The same as with four methods: the Presence Pearson r of these pairs (0.42–0.78) is below its 0.85 call, so it supports none of them in RNA.
  - DNA: P01–EE05 from step 2, contig_0000623–EE07 from step 3 (the Presence Pearson supports it, r 0.90 over 20 libraries), the others at no step.
- **Single-method line style:** with exactly one method switched on, line width and opacity follow its score in five classes (`network_site.SCORE`): FlashWeave |weight| < 0.3 / 0.3 / 0.4 / 0.5 / ≥ 0.7; SPRING and SPIEC-EASI stability < 0.4 / 0.4 / 0.6 / 0.8 / ≥ 0.9; Spearman |ρ| < 0.3 / 0.3 / 0.4 / 0.5 / ≥ 0.6; Presence Pearson |r| < 0.87 / 0.87 / 0.9 / 0.93 / ≥ 0.96 (its edges start at 0.85). A legend under the colour key shows the classes; with several methods it shows the width per number of supporting methods.

## Rebuild of 2026-11-09 (30/20 node rule + pool expression rule)
- Nodes: DNA 551 (454 prokaryotes, 25 EE, 59 vOTUs, 13 plastid), RNA 602 (496 / 25 / 63 / 18). Prokaryote flags (RNA layer): 54 fragment, 67 contamination 10–20%, 5 both.
- Run time of `run_all.sh` (16-core machine, 602 RNA nodes): nodes.py 25 s, SPRING both layers 13 min, SPIEC-EASI 30 s, FlashWeave 44 s, the rest under 1 min; 15 min in all.
- Links: FlashWeave 888 / 1,870; SPRING 7,308 / 10,356; SPIEC-EASI 4,114 / 6,142; residual Spearman 5,477 / 82,897; Presence Pearson 6,842 / 2,668 (4.5% / 1.5% of pairs; SPRING 7,308 / 10,356, so still at or below SPRING); consensus 2,065 / 9,422; sign conflicts 68 / 306 (DNA / RNA).
- Possible same genome: 85 pairs (Pelagibacter, UBA4592, Sulfurimonas, ...); all 85 are RNA consensus links, 71 are DNA links (45 in the DNA consensus). The four Cycloclasticus MAGs are not among them.
- Steps (DNA / RNA links): 70 / 463, 236 / 1,076, 929 / 3,535, 1,046 / 3,846, 2,065 / 9,422 (consensus), 3,798 / 9,987, 5,306 / 14,511, 19,109 / 40,216.
- Known pairs: RNA P01–EE05 from step 1; contig_0000362–EE04 and contig_0000564–EE07 from step 5; contig_0000610–EE22 from step 6; contig_0000623–EE07 from step 7; the virophage at no step. DNA: P01–EE05 and contig_0000623–EE07 from step 3.
- Page: positive / negative / both link filter (also applied to the counts and the known-pair table, where ○ marks a pair found as a negative link); "find genome" box (id or any part of the name; rings and selects matches, scrolls to the map); Cycloclasticus MAGs always labelled (`network_site.FEATURED`); sample picker (dropdown, ◀ ▶, arrow keys, ▶▶ one sample per second) sizing points by the genome's level in that library from `site/network_{layer}_samples.json` (RNA reads per kb per million of the robust profile, DNA parts per million of aligned bases; 2 significant digits; RNA 472 kB, DNA 87 kB, loaded on first use). Order: winter > control, diesel > seawater, sea ice top, other ice, sea ice bottom > day; non-pool libraries last.
- Sizes: network.json 92 kB, network_dna.json 446 kB, network_rna.json 807 kB.

## Illumina DNA layer (2026-11-09, layer key `mgdna`)
- **Data:** danaSeq illumina_rna in DNA mode (`--rrna_method none`, human removal, strobealign, ≥ 80% aligned) of the CMO Illumina metagenomes against the same coassembly (Grex `~/scratch/project_rnaseq/marine_mg` → `data_mg/`). Loader `common.illumina_dna_samples()` / `illumina_dna_contig_matrix()`: lanes (305) summed per library (169); the 7 QC-flagged libraries (`excluded('metagenome')`) are dropped, which already removes the four failed ones (< 0.2 M raw reads); the ≥ 1 M pairs-on-contigs floor drops nothing more. 162 libraries remain; the layer uses the 85 pool libraries (CTL/EXP, seawater or sea ice, with a day; incubations, MOS cultures, sites, foam and enrichments left out, as for nanopore). Each library is its own sample (replicate libraries are not summed). Median 72% of non-human reads mapped (pool libraries 78%).
- **Metric:** read pairs on the node's contigs; library size = all contig-assigned pairs; the methods get the read pairs as counts (`unit` 1). With fixed-length 2 × 150 reads this is the nanopore metric (depth × length share), so both DNA layers use one metric, summed over contigs (DNA is not spiky like RNA; `euks/mags/scripts/mgdna.py`). Sample picker: per million mapped read pairs.
- **Nodes:** same candidates and pool-expression rule; prevalence ≥ 25% and ≥ 1 per million (EE, giant viruses, virophage, organelle groups ≥ 10%). Illumina DNA is deep, so almost every candidate passes: 598 nodes (496 prokaryotes, 25 EE, 61 vOTUs, 16 organelle).
- **Results** (links per method; consensus = ≥ 2 methods, same sign, one design-aware):

| Layer | Nodes | FlashWeave | SPRING | SPIEC-EASI | Residual Spearman | Presence Pearson | Consensus | Sign conflicts |
|---|---|---|---|---|---|---|---|---|
| Nanopore DNA (44) | 551 | 890 | 7,308 (λ 0.042) | 4,114 (λ 0.144) | 5,477 | 6,842 | 2,066 | 68 |
| Illumina DNA (85) | 598 | 1,474 | 9,730 (λ 0.026) | 6,677 (λ 0.054) | 89,047 | 6,942 | 10,692 | 736 |
| RNA (128) | 602 | 1,871 | 10,356 (λ 0.033) | 6,142 (λ 0.088) | 82,897 | 2,668 | 9,421 | 307 |

  - Illumina DNA steps 1–8: 607, 1,294, 3,758, 4,743, 10,692 (consensus), 11,608, 25,461, 60,330 links. 34,655 Spearman-only links below |ρ| 0.4 are left out of `site/network_mgdna.json` (1.1 MB, the largest network file; samples file 204 kB).
  - The nanopore DNA and RNA layers were rerun with everything else (`run_all.sh`); only FlashWeave moved slightly (DNA 888 → 890, RNA 1,870 → 1,871 links; its convergence heuristic is not deterministic), so the consensus moved by one link each (2,065 → 2,066, 9,422 → 9,421). Known pairs unchanged.
  - Possible same genome: all 85 pairs are Illumina DNA links, 77 in its consensus.
- **Known pairs (Illumina DNA):** P01–EE05 from step 1 (all five methods); giant virus contig_0000564–EE07 from step 5 (consensus: residual Spearman ρ 0.39 + Presence Pearson r 0.91 over 57 libraries). Not found at any step: contig_0000623–EE07 (Presence Pearson r 0.82, below its 0.85 call; residual Spearman ρ 0.27), the virophage contig_0049668–EE07 (r 0.54 over 58 libraries; Spearman ρ 0.29, q 0.026, below the 0.4 export floor), contig_0000362–EE04 (ρ 0.29) and contig_0000610–EE22 (Spearman negative, ρ −0.31). So in DNA the virophage does not follow EE07 closely; in RNA it is not found either.
