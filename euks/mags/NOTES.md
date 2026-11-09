# Eukaryotic MAGs on marine-coassembly-flye (2026-11-08)

## Scripts

Run in this order: `survey.py` → `compact.py` → `final.py` → `expr.py`.

- **`survey.py`**: per-bin stats for the user's listed bins; overlap matrices (Jaccard and containment); first grouping by containment ≥ 0.5.
- **`compact.py`**: t-SNE compactness against prokaryotic MAGs, fragmentation, footprints, and the density-gap refinement.
  - Spatial clusters come from DBSCAN(eps 1.5, min_samples 20, weighted by contig kb).
  - Each cluster holding ≥ 5% and ≥ 2 Mbp of a group becomes its own genome.
  - Smaller clusters join a genome within 24 units (10% of the map width) and ≤ 3 GC points; the rest are trimmed.
- **`final.py`**: merges groups that sit in each other's footprint (≥ 50% of the smaller one inside, GC within 2 points). Also finds uncaptured Tiara-eukarya blobs and gathers taxonomy hints.
- **`expr.py`**: reads per group from full-depth contig counts (any reads on the contig, so introns are fine) and CDS counts.

## Main outputs (`out/`)

- `euk_groups_final.tsv`
- `contig_euk_group.tsv` (contig → E##)
- `uncaptured_blobs*.tsv`
- `group_expr_*.tsv`
- PNGs: `tsne_final.png`, `tsne_euk_region.png`, `tsne_EUK*.png`

## Key facts

- The listed bins form 20 groups.
  - Core genomes (purity ≥ 0.9, ≥ 4 Mbp): E01–E07, E10, E12, E13.
  - Purity about 0.5, probably mixed or not eukaryotic: E08 (`metabat_4578`/`vamb_3866` = `semibin_630`), E09 (`comebin_1009`; 15% Tiara archaea, Whokaryote 98% prokaryote), E11 (`vamb_1186`, depth 0), E20 (`binette_1014`).
- `semibin_1301` is three density-separated genomes with distinct GC: E02 (49.3%), E06 (52.7%) and E13 (53.9%).
- `semibin_1280` holds two genomes: E01 (46.9% GC) and E10 (37.3% GC; the `comebin_1991` / `vamb_3134` / `metabat_10781` / `binette_996` family).
- Missed eukaryotic content: Tiara-eukarya blobs outside the listed bins.

| Blob | Size | Main bin | Note |
|---|---|---|---|
| U5 | 52 Mbp | `semibin_1338` | Same spot as E05 (`comebin_3937`) |
| U4 | 17 Mbp | `semibin_1317` | 59% GC |
| U3 | 14 Mbp | `semibin_1315` | |
| U0.4 | 9.6 Mbp | `semibin_1384` part | Whokaryote 0.99 |

  - About 830 Mbp of diffuse eukaryotic sequence sits in the catch-all bins `semibin_1397`/`1398`/`1384`.

## Size limits in danaSeq `mag_analysis`

- `nextflow.config:109` (`checkm2_max_bin_bp`) and `:122` (`gtdbtk_max_bin_bp`) are both 50 Mbp. They are applied at `modules/binning.nf:997` and `modules/phylogeny.nf:37`. These only skip CheckM2 and GTDB-Tk on bins above 50 Mbp.
- `nextflow.config:124`: `binette_max_bin_bp = 10 Mbp`, applied at `modules/binning.nf:640`. Bins larger than this never reach Binette, which is why the Binette eukaryotic bins are 2–7 Mbp fragments.
- DAS Tool scores bins with prokaryotic single-copy genes, so it rarely keeps eukaryotic bins. Only `comebin_3937` and `comebin_4527` got through.

## v2 (2026-11-08): added candidates, stable IDs, wiring tables

**Candidates.** The bin list is now the config `LIST` at the top of `scripts/survey.py`. Added: `semibin_1338`, `semibin_1315`, `semibin_1317`, `semibin_1384`.

**Rules added**
- **Genome seeds** (`compact.py`): a seed must also be ≥ 50% Tiara-eukaryotic by length. This drops the bacterial parts of `semibin_1384`, `vamb_1186` and `binette_1014`.
- **Stable IDs** (`final.py`, reference `out/v1/`):
  - A new group keeps an old E## if it holds ≥ 50% of that group's bp and the old group is ≥ 1/3 of the new one.
  - Smaller old groups are "absorbed" and noted.
  - New groups take the next free number.
- **`network_ok`** (`finalize.py`): all of
  - Tiara purity ≥ 0.9;
  - ≥ 4 Mbp;
  - median contig-read share ≥ 10 per million;
  - t-SNE RMS radius ≤ 12.
- **Status**: one of network / low RNA / <4 Mbp / doubtful.

**Changes from v1**
- E05 grew from 27.9 to 75.1 Mbp by merging with `semibin_1338`.
- New groups:
  - E21: 145 Mbp from the diffuse eukaryotic part of `semibin_1384`. It absorbed E19 (the `metabat_4788` fragment). Spread 20 units, so doubtful.
  - E22: `semibin_1317`, 59% GC. Purity 0.89, so doubtful.
  - E23: `semibin_1315`. Network-ready.
  - E24: a 0.27 Mbp fragment. Doubtful.
- Dropped by the purity rule: E11 (`vamb_1186`), E19 (absorbed into E21), E20 (`binette_1014`).
- All other groups unchanged.

**Outputs for wiring**
- `out/contig_euk_group.tsv`: contig, euk_group, source_bins, tiara, gc, length.
- `out/euk_groups.tsv`.
- `out/euk_groups_tsne.png`.

## Wired into the site (2026-11-08)
- **Pipeline:** `scripts/bin_domains.py` copies the wiring tables into `data/`.
  - `data/euk_groups.tsv`: one row per group.
  - `data/contig_euk_group.tsv.gz`: contig → E## plus its Tiara call.
  - It also adds `bin_domain.tsv:euk_group`: the group holding ≥ 50% of a DAS Tool eukaryote bin's bp. Result: `comebin_3937` → E05, `comebin_4527` → E03. `lorbin_494_sub` is in no group and stays "Eukaryote (Tiara)".
- **One group per contig:** 60 contigs had been assigned to two overlapping groups (mostly E21 with E17 or E01). `contig_euk_group.tsv` is now deduplicated upstream: each goes to the larger group, so E12, E15, E17 and E24 shrank slightly. `bin_domains.py` keeps the same rule as a guard.
- **Helpers in `scripts/common.py`:**
  - `euk_group(contig)` and `euk_group_info(gid)`;
  - `is_euk_group(x)`;
  - `euk_group_label(gid)`: "Eukaryote E02 (Ochromonadales)", "Eukaryote E05", or "Eukaryote E21 (mixed group)" for doubtful groups;
  - `contig_host_label(contig, fallback)`: for unbinned contigs. Contigs in a doubtful group take its label only if Tiara calls the contig itself eukaryotic.
  - `special_label` maps eukaryote bins to their group's label.
- **E03** has two dense lobes about 10 t-SNE units apart. It passes the compactness rule (RMS 5.4) but may be two related genomes. It is kept as one group for now.
- **Reproducibility check:** rerunning survey → compact → final → expr → finalize on 2026-11-08 reproduced all 21 groups exactly: the same bp, contigs, purity, status and expression.
- 2026-11-08: the figure has a unit grid (t-SNE 1/2) and a third panel. That panel shows the eukaryotic region with contigs ≥ 2 kb sized by cumulative expression (per-library contig reads per million, summed over the 128 analysed libraries); labels give group totals.
- Status "low RNA" was renamed "little CMO RNA": the assembly is about 50% non-CMO samples, so this is not a quality flag.
- Contigs claimed by two groups now go to the larger group. 60 contigs were affected.

## Blob capture, viruses, differential abundance (2026-11-08)

**Run order:** survey → compact → final → expr → finalize → kernel → da → register_blobs → plot.

**`kernel.py`** captures blobs at the user's seed points (config `SEEDS`) with an expression/coverage kernel.
- **Profiles:** RNA = per-library contig reads per million (128 libraries); DNA = per-sample depth × length share from Grex `mag_marine/mapping/depths.txt` (72 CMO nanopore libraries; 30 of the 102 are not in the table).
- **Score:** Gaussian spatial kernel (σ 3) × mean Pearson r with the blob's log profile. Members score ≥ 0.5; the profile is recomputed until membership is stable.
- **Why not density alone:** pure density clustering (DBSCAN, HDBSCAN) either merges the whole eukaryotic region or splits blobs (see `seeds_hdbscan.png`).

**`da.py`** runs per-blob DA in RNA and DNA.
- Contrasts: diesel vs control after dosing, per season; BACI 2025-26 seawater; ice vs water, paired.
- BH correction within each contrast × data type.
- DNA has no 2025-26 diesel-vs-control or BACI test: the control pool's nanopore libraries stop at day 3.

**`register_blobs.py`** is the rule for turning blobs into new E##/V## groups.
- Diffuse = surrounding Tiara-eukarya contigs outweigh the blob AND follow its profile (median r ≥ 0.6).
- Results:
  - New eukaryotic groups E25–E30.
  - V02 (new virus, (−48, −28)).
  - B01 = E22; B10 = V01.
  - B04 is a fragment.
  - B02 and B11 are prokaryotic (B11 = `comebin_277`).

**Viruses:** groups with ≥ 50% of bp geNomad 'Viruses' get V## IDs (`virus_ids.tsv`) and are kept out of the eukaryotic tables. E09 became V01.

## Split rule and shared ID registry (2026-11-08)

**`scripts/eukshared.py`** holds two shared pieces.
- **Split rule (`split_parts`):** a group splits into the two parts of a 2-component GMM on contig t-SNE coordinates when all of these hold:
  - Ashman's D > 3;
  - the centres are > 5 units apart;
  - the minor part is ≥ 2 Mbp and ≥ 10% of the group.
- **ID allocator:** one allocator, `out/id_registry.tsv`, keyed by `split:<parent>:<k>`, `blob:<seed>` and `group:<refined group>`.
  - Used by `final.py` (splits and new groups) and `register_blobs.py` (blobs; the split test also runs there).
  - IDs are stable across reruns.

**Effect on current groups:** only E03 split.
- E03 (24.8 Mbp, the larger lobe) and E31 (7.6 Mbp, "split from E03", centres 10.3 units apart). Both lobes are 100% `comebin_4527`.
- Kernel blobs have D 1.9–3.7 but are compact (centres < 5 units apart), so none split.
- E21 has D 2.7 and stays whole.

## Site wiring after the split, blob and virus update (2026-11-08)
- **New data tables:** `scripts/bin_domains.py` now also writes `data/virus_groups.tsv` and `data/contig_virus_group.tsv.gz`.
- **New helpers in `common.py`:**
  - `virus_group`, `virus_group_label` ("Virus V02"), `is_virus_group`;
  - `group_expression()`: E## and V## × libraries from `group_expr_samples.tsv`. IDs are renamed via `virus_ids.tsv`, and kernel-blob groups use the blob's RNA abundance through `blob_registry.tsv`;
  - `group_da()`: per-group DA from `da_blobs.tsv` through `blob_registry.tsv`.
  - `contig_host_label` labels V## contigs "Virus V##" before checking eukaryotic groups.
- **`comebin_4527` → E03:** it holds 77% of the bin's bp.
- **V02 is probably plastid DNA, not a virus:**
  - 98.7% of its 4.95 Mbp is Tiara organelle;
  - its contigs sit in the plastid MAGs `metabat_1683`, `lorbin_892`, `metabat_11799`, `lorbin_663` and others;
  - Kaiju and 16S say Cyanobacteriota;
  - it carries psbA/D, psaA/B, petB/D and rbcL: on the Light page "Virus V02" hosts 15–30% of those transcripts.
  - geNomad calls 78% of it viral ("Heunggongvirae"). The V## rule (≥ 50% geNomad viral) should probably exclude groups that are ≥ 50% Tiara organelle.
- **Viruses page:** shows an organelle column and a caution for V02.

## STATE FOR RESUMING (2026-11-08, end of session)

### Pipeline order (all in `scripts/`, Python `/data/miniforge3/bin/python3`)
1. `survey.py` (candidate config `LIST`) → `compact.py` → `final.py` (split rule + ID registry) → `expr.py` → `finalize.py` (virus rule, statuses).
2. `kernel.py` (seed config `SEEDS`) → `da.py` → `register_blobs.py` (organelle → P/M, virus → V, eukaryote → E) → `plot.py`.
3. `hull_panels.py` → **`hulls.tsv`** (curated, hand-picked; never regenerated) → `hull_mags.py` → **`merge_mags.py`** (reruns `hull_mags.py` when the merges change) → `hull_expr.py` → `plot_hulls.py` → `heatmap.py`.
4. BUSCO and identity: `busco_table.py` (runs BUSCO on any MAG whose membership has no cached run) → `busco_dups.py` → `identify.py` → `coding.py` (adds BUSCO, coding and `genome_status` to `ee_taxonomy.tsv`). Then the site: `scripts/bin_domains.py` and the page scripts.
5. `organelle.py`: exploratory, see Pending. Do NOT run it as part of the pipeline: it rewrites the organelle tables with a flawed split.
- Always run `final.py` before `finalize.py`. `finalize.py` rewrites `contig_euk_group.tsv` from bin groups, and `register_blobs.py` appends the kernel blobs afterwards.

### Rules (see sections above for details)
- **Eukaryote groups (E##):**
  - Built by containment + density gaps (DBSCAN eps 1.5).
  - Seeds must be ≥ 50% Tiara eukarya.
  - GC/24-unit merge; footprint merge.
  - Split rule (`eukshared.split_parts`): Ashman D > 3, > 5 units apart, minor part ≥ 2 Mbp and ≥ 10%, **and** the parts' RNA profiles differ (Spearman ρ < 0.9 between the parts' summed per-library RNA; guard added 2026-11-08).
- **Viruses (V##):** ≥ 50% geNomad viral AND < 50% Tiara organelle.
- **Organelles (P##/M##):** ≥ 50% Tiara organelle. P vs M by marker genes, because Tiara's second-stage call is not in the contig table.
- **IDs:** shared registry `out/id_registry.tsv`.
  - V02 is retired: it was a plastid, now P01. The keys `retired:V02` and `retired:VV02` keep those numbers from ever being reused.
- **Hull MAGs:**
  - EE##/VV## = point-in-polygon on `hulls.tsv`, contigs ≥ 1 kb.
  - EE hulls exclude Tiara bacteria/archaea.
  - Overlaps go to the source group, else the nearer centroid.
  - Single-hull sources map E05 → EE05; extra E21 cores → EE32–EE36.
  - `network_ok` (EE) = purity ≥ 0.9, ≥ 4 Mbp, median ≥ 10 reads per million.

### Output files (`out/`) and their state
- **Current, wired by the coordinator:** `contig_euk_group.tsv` and `euk_groups.tsv` (E01–E08, E10, E12–E18, E21–E31). Contigs are unique.
- **Current:**
  - `virus_groups.tsv` / `contig_virus_group.tsv` (V01 only).
  - `organelle_groups.tsv` / `contig_organelle_group.tsv` (P01, unsplit, 335 contigs, 4.95 Mbp).
  - `hull_mags.tsv` / `contig_hull_mag.tsv` / `hull_mags_da.tsv` (EE01–EE36 + VV01; 85,878 contigs, unique).
  - `euk_groups_tsne.png`, `euk_hull_mags_tsne.png`.
  - `hulls/panels/*.png` (drawing panels), `hulls/verify/*.png`, `hulls/verify_expr_sheet_*.png`.
- **Snapshots:** `*.before_split.tsv` (pre-E03-split eukaryote table; pre-split organelle tables = the current ones). `out/v1/` is the ID-continuity reference for `final.py`.
- **Marker table:** `organelle_markers.tsv`, the plastid/mito marker copies on P01 by bin.

### Pending (not done; user deferred)
1. **Plastid P01 split.**
   - P01 holds about 21 plastid genomes of markers (rbcL 33, psaA 29, psaB 34, psbD 23, petB 45, atpB 8; mito cox1 5, cob 1). Most copies sit on unbinned contigs (rbcL 23, psaB 25).
   - Nine organelle bins inside it are single plastids of 100–130 kb, 1 marker copy each: metabat_1343, lorbin_603, metabat_11799, lorbin_507, lorbin_663, lorbin_613, lorbin_732_sub, lorbin_892, metabat_7332_sub (16S host Thalassiosira pseudonana).
   - Plan:
     - one P per organelle bin;
     - split the unbinned remainder by RNA+DNA profile co-variation (average-linkage, r ≥ 0.85; keep clusters ≥ 50 kb carrying rbcL/psaA/psaB);
     - the parent ID goes to the largest part.
   - The `organelle.py` draft gives the parent ID to the largest bin and lumps the remainder into it. Fix that before use.
2. **Host links.** Spearman of each plastid's RNA profile against the E## groups. Draft top matches: E05 for most bins (ρ 0.8–0.97), E07 for lorbin_663, E21 for metabat_1343/metabat_11799. Lineage agreement is not checked yet.
3. **Coassembly-wide organelle check.** Tiara-organelle contigs ≥ 1 kb: 313 in P01 (4.9 Mbp), 36 in 16 other organelle bins (1.6 Mbp; plastid MAGs elsewhere on the map, e.g. metabat_9288_sub, lorbin_635, lorbin_511, metabat_7265, metabat_8525), 105 in non-organelle bins (1.2 Mbp), 9,808 unbinned (53.7 Mbp, location not yet mapped).
4. **Hulls to re-check with expression (user asked whether hulls used expression; they did not).** Hulls were drawn on spatial density with contigs sized by length. Review the expression check sheets `hulls/verify_expr_sheet_*.png` for:
   - EE01 (median expression 1 → 104 per million after gaining 2.3 Mbp: probably an expressed neighbour inside the hull);
   - EE12 (+6.1 Mbp), EE14 (+6.0), EE30 (+14.6), EE31 (+9.5), EE27 (+7.6): large gains;
   - EE36 (E21F) and EE13: small or uncertain;
   - EE08 (6.5 Mbp of prokaryotic contigs excluded).
5. **M## series:** no mitochondrial group exists yet (the rule is ready).

## Site wiring of the hull MAGs (2026-11-08, PR #17)
- **EE## (and P##) are the site's eukaryotic and organelle labels.** E##/V## stay in `data/` as provenance only; each EE records its `source_group`.
- **Expression:** `scripts/hull_expr.py` (run after `hull_mags.py`) writes `out/hull_mags_expr_samples.tsv`: EE/VV/P × 128 libraries, reads per million contig-assigned reads. It uses the same metric as `expr.py`, and its medians equal `hull_mags.tsv:expr_ppm_median` exactly. DA comes from `out/hull_mags_da.tsv`.
- **`scripts/bin_domains.py`** writes:
  - `data/hull_mags.tsv`, with `status` (same thresholds as the E## groups) and a rule-based `caution`:
    - "uncertain: hull may include an expressed neighbour": median expression ≥ 10× the source group's and ≥ 10 per million. Only EE01 (1.1 → 104).
    - "hull added X Mbp to the Y Mbp source group": gained ≥ 4 Mbp and ≥ 50% of the source. EE12, EE14, EE16, EE26, EE27, EE29, EE30, EE31.
    - "hull drawn on a sparse or small core": hull note starts "UNSURE". EE12, EE13, EE36.
  - `data/contig_hull_mag.tsv.gz`, `data/organelle_groups.tsv`, `data/contig_organelle_group.tsv.gz`;
  - `bin_domain.tsv:hull_mag` (comebin_3937 → EE05, comebin_4527 → EE03, lorbin_494_sub → EE30).
- **Viral genome-group rule (VV## and V##):** a viral group needs ≥ 50% of its bp in geNomad Nucleocytoviricota or another giant-virus / eukaryotic-virus class (Megaviricetes, Maveriviricetes, Pokkesviricetes, Polintoviricetes, Alsuviricetes, Pisoniviricetes, Revtraviricetes).
  - Tailed-phage (Caudoviricetes) clusters are many phages that cluster together on composition, not genome groups.
  - V01/VV01 has 2,165 contigs and 20.9 Mbp, is 97% Caudoviricetes, spans 1–405× depth, and its largest contig is 569 kb. It fails the rule (0% giant-virus classes), so **no VV## exist**.
  - `retired:VV01` is in `out/id_registry.tsv`; nothing was renumbered. Its contigs keep the per-contig vOTU handling on the Viruses page.
  - `hull_mags.py` now applies the rule. The current `out/hull_mags.tsv`/`contig_hull_mag.tsv`/`hull_mags_da.tsv` predate it and still list VV01; `bin_domains.py` drops it when wiring `data/`, and `plot_hulls.py` leaves retired ids out of the figure.
- **`contig_host_label` order:** VV## (none at present) → P## ("Plastid P01", counted as plastid like `is_plastid`) → EE## (doubtful EE only if Tiara calls the contig eukaryotic) → fallback.

## EE30 re-check with expression (2026-11-08)
- **Why:** the user saw EE30 sitting inside a larger dense cloud of expressed contigs (bottom panel of `euk_hull_mags_tsne.png`).
- **Tool:** `scripts/hull_check.py <hull> x0 x1 y0 y1` draws a gridded zoom (contigs ≥ 2 kb sized by cumulative RNA, coloured by Tiara) and a second panel coloured by each contig's mean RNA/DNA Pearson r with the hull's reference profile. The first run saves the current members as the reference (`out/hulls/check/<hull>_ref.tsv`), so a later polygon edit is judged against the old core. It also prints median r by contig length for reference, new and outside contigs.
- **Finding:** in rings around the EE30 centre (−1.9, −38.6), r at matched contig length is the same out to radius ~4.5 and drops past 5. GC is flat (42.9–43.0%). The raw mean r falls earlier only because the halo contigs are shorter. The anticorrelated streak at (4, −42) and E31 below are kept out.
- **Polygon:** EE30 in `hulls.tsv` is now a 13-vertex polygon of radius ~5. The old disc was (−4.5,−38.5) (−3.5,−36.3) (−1.5,−35.8) (0,−36.5) (0.6,−38.5) (0,−40.5) (−1.5,−41.2) (−3.5,−40.5).

| EE30 | before | after |
|---|---|---|
| Size | 36.95 Mbp, 2,721 contigs | 71.77 Mbp, 7,855 contigs |
| Tiara eukarya | 0.998 | 0.998 |
| Median expression | 342 per million | 684 per million |
| t-SNE RMS | 1.64 | 2.84 |
| N50 | 27.5 kb | 13.3 kb |

- **Added contigs:** 34.8 Mbp in 5,134 contigs, none of them in another MAG before.
  - Summed profile against the old core: RNA ρ 0.96 and DNA ρ 0.80 (neighbouring MAGs reach ρ 0.61–0.85).
  - Contig r with the old core profile at matched length: 0.25/0.34/0.43/0.62 (1–5, 5–10, 10–30, > 30 kb), against 0.26/0.34/0.45/0.71 for the core and 0.04/0.12/0.23/0.26 for contigs outside.
  - Conclusion: the added contigs are the same genome.
- **Proteins:** the expanded MAG is Chaetoceros (79% of MarFERReT votes).
- **Caution:** EE30 now gets the "hull added 48.4 Mbp" caution. That rule fires on size alone.
- **Downstream reruns:** `hull_mags.py`, `hull_expr.py`, `plot_hulls.py` and `heatmap.py`. `hull_mags.py` now also drops the retired VV01 from `out/`.

## Taxonomy of the EE MAGs (2026-11-08, `scripts/identify.py`)
- **Evidence:** existing danaSeq outputs, fetched once from `mag_marine` with `scripts/fetch_grex.sh` into `out/identify/grex/` (not committed):
  - MetaEuk proteins against MarFERReT (`eukaryotic/marferret/marferret_proteins.tsv`). The pipeline left the taxonomy empty, so hits are mapped to NCBI lineages locally (MarFERReT taxonomies plus the BLAST taxdump).
  - barrnap + vsearch/SILVA rRNA genes (`taxonomy/rrna/rrna_genes.tsv`).
  - Kaiju from `contigs_full.pkl`, as a cross-check only. Kaiju classifies 1–55% of EE bp and mostly says Pseudomonadota, so it is uninformative for these MAGs.
- **Rules** (the docstring has the details):
  - The protein call is the deepest NCBI rank with ≥ 50% of alignment-length-weighted votes and ≥ 20 proteins. When the median identity is < 50%, the call goes no deeper than class.
  - rRNA counts only with ≥ 2 nuclear rRNA contigs. Single rRNA contigs are often misbinned, e.g. a fungal 28S in EE30.
  - Confidence:
    - high: support ≥ 0.6, plus rRNA agreement or median identity ≥ 90%;
    - medium: support ≥ 0.5 and not contradicted; or protein ≥ 0.3 corroborated by rRNA; or rRNA-only from ≥ 3 contigs;
    - low: a best guess only.
- **Output:** `out/ee_taxonomy.tsv` (with an empty `completeness` column; BUSCO is being run separately into `out/busco/`), copied to `data/ee_taxonomy.tsv`.
- **Site wiring:** `hull_mags.tsv` gains `taxonomy`, `taxonomy_rank` and `taxonomy_confidence`. `hull_mag_label` names a MAG only at high or medium confidence. The old 18S `taxonomy_hint` is no longer used in labels.
- **Calls:**
  - Diatoms: EE05 Minidiscus (Thalassiosirales; its linked plastid's 16S says T. pseudonana), EE30 Chaetoceros, EE03 Attheya, EE29/EE31 Chaetocerotales, EE27 Thalassiosiraceae, EE26 Staurosira, and EE21 and EE32–EE35 pennates.
  - Green algae: EE22 Micromonas, EE23 Bathycoccus, EE25 Picochlorum.
  - EE02 and EE28 Chrysophyceae; EE04 Cercozoa (rRNA); EE07 Acanthoecida (choanoflagellate).
  - EE12, EE14 and EE17: Calanidae (Calanus copepods), all "little CMO RNA". EE10: Metazoa (94% of kingdom votes, no close reference).
  - EE01: Laminariales by proteins (kelp-like). Its two 18S/28S contigs say Cercozoa, but they sit on the hull's west edge.
- **EE08:** a ciliate (Ciliophora 0.82 of protein votes, Oligohymenophorea 0.64, median identity 43%, so distant from references). It is not a metazoan: Metazoa gets < 4% of votes.
  - 95% of its RNA sits on 5 contigs (24 kb) present in 123–128 of 128 libraries.
  - Those contigs carry tubulin, ubiquitin and HSP90 with 72–96% identity to many lineages. They are conserved, highly expressed genes, probably taking cross-mapped reads from many eukaryotes.
  - The rest of the genome has a median of 0.8 reads per million, so EE08's heatmap RNA row is not the genome's own expression.

## Split guard, merge rule and relabel rules (2026-11-08, user-approved)
- **Split guard** (`eukshared.split_parts`): a density split now also needs the parts' RNA profiles to differ, Spearman ρ < 0.9 between the parts' summed per-library RNA (reads per million contig-assigned reads, 128 libraries; the `hull_expr.py` metric, `eukshared.rna_profile`).
  - E03's two lobes have ρ 0.996, so E03 is whole again. `final.py` retires split IDs a rerun no longer produces: `retired:E31 (split of E03 undone by the RNA split guard, 2026-11-08)`.
  - The E31 hull keeps its id EE31: `hull_mags.py` resolves a hull whose source group was retired by the guard to the parent (`eukshared.resolve_source`: E31 → E03).
  - Reran final → expr → finalize → da → register_blobs → plot; only the 705 E31 contigs changed (to E03). Kernel blobs, viruses and organelles are unchanged.
- **Merge rule** (`scripts/merge_mags.py`, after `hull_mags.py`). A pair of EE MAGs merges only if all hold: RNA ρ ≥ 0.95; DNA ρ ≥ 0.8 (`eukshared.dna_profile`, CMO nanopore depth × length share); skani AF < 10% both ways (`skani dist --min-af 0`; no output = below the 80% ANI screen = AF 0); BUSCO on the union D < 5% and C above the best part.
  - Transitive in decreasing RNA ρ; a growing group is re-tested as a whole (summed profiles, skani of group FASTAs, BUSCO of the union).
  - The lowest EE number is kept; others get `retired:EEnn (merged into EEmm, <date>)`. `hull_mags.py` applies these keys (`eukshared.merged_into`); merged MAGs list all hulls/source groups (`E03+E29`) and a `merged_from` column. `merge_mags.py` reruns `hull_mags.py` when the merge set changes.
  - Every pair goes to `out/mag_merge_tests.tsv`. skani runs for pairs with RNA ρ ≥ 0.8 or the same source group; union BUSCO runs for candidates and near-misses (same source group, or RNA ρ ≥ 0.9 failing one rule).
- **BUSCO cache** (`scripts/buscolib.py`): runs are keyed by the md5 of the sorted contig ids (`out/busco/runs/<label>_<md5>/`, with `members.txt`); older runs under `out/busco/run_*` are adopted when their input FASTA still exists, predates the summary and holds exactly the MAG's contigs. FASTAs are cut from `ee_contigs.fa.gz` + `ee_contigs_add1.fa.gz`.
  - `busco_table.py` writes `out/busco/busco_eukaryota_odb12.tsv` for the current MAG set (newest matching run per MAG; runs BUSCO when none matches). `busco_dups.py` reads the chosen runs.
- **Relabel rules** (`coding.py` → `ee_taxonomy.tsv:genome_status`; `bin_domains.py` → `data/hull_mags.tsv`; `common.hull_mag_label`):
  - `doubtful, not a genome`: BUSCO C < 1% AND coding < 2% AND not a metazoan call (protein kingdom vote Metazoa ≥ 0.5, or a metazoan taxon in the call). Status doubtful, not network-ready, label "EE01 (doubtful)".
  - `multi-genome (≥ 2 species)`: BUSCO D ≥ 10% AND median duplicate-copy identity < 90%; `(≥ 2 strains)` for 90–98%. Stays network-ready with a caution; label e.g. "EE30 (Chaetoceros spp., ≥ 2 species)".
  - Hulls whose note in `hulls.tsv` contains `RECHECKED <date>` (E30: expression re-check; E12: copepod lineage panel) are exempt from the "hull added X Mbp" caution.
- **Heatmap:** a MAG whose 5 most-read contigs hold ≥ 90% of its RNA reads is labelled `EE##*` and its RNA row is hatched (`out/hull_mags_rna_flags.tsv`).
- **Result (2026-11-08):** merged EE29 → EE03 (RNA ρ 0.979, DNA 0.854, AF 1.6/4.4%, BUSCO 55.8 [2.3] + 5.4 [0] → 59.7 [2.3]) and EE13 → EE06 (ρ 0.950, DNA 0.827, AF 0.2/2.9%, 12.4 + 0.8 → 13.2 [0]; borderline on RNA and one extra BUSCO). Not merged: EE03 + EE31 (RNA ρ 0.53; the union would give 64.3 [3.9]). EE31's summed RNA is 83% on its 5 most-read contigs (two 10-kb contigs, contig_0221272 and contig_0193859, in the 9.5 Mbp the hull added); without its top 5 contigs ρ with EE03 is 0.97, and its E31 core alone has ρ 0.99. E21 pieces: best RNA ρ 0.82 (EE33/EE35), none pass. All pairs: `out/mag_merge_tests.tsv`.
- **RNA flags:** EE01 and EE04 have 83% and 94% of their reads on one 18S/28S rRNA contig each; EE08 on conserved-gene contigs; EE10 on two contigs. Their RNA medians (e.g. EE04 3,107 per million) are not genome expression.
- **Distant genus calls read "-like"** (user rule 2026-11-08; `common.taxon_display`, written by `scripts/bin_domains.py` as `data/hull_mags.tsv:taxonomy_display` and used by `hull_mag_label`, the Eukaryotes table and `lineage.json`): a genus call with median protein identity (`protein_pid_median`) < 90% shows the deepest of order/class with protein support ≥ 0.5, plus the genus as "-like": EE25 "Trebouxiophyceae (Picochlorum-like)", EE03 "Chaetocerotales (Attheya-like)", EE26 "Fragilariales (Staurosira-like)", EE30 "Chaetocerotales (Chaetoceros-like)". Genus calls at ≥ 90% keep the plain genus (Micromonas, Bathycoccus, Minidiscus). Inside a label the brackets become commas: "Eukaryote EE30 (Chaetocerotales, Chaetoceros-like, ≥ 2 species)".
- **Published Hudson Bay context:** `euks/literature_context.tsv` (curated by hand; `#ref` lines hold citations) is read by `euks/scripts/euks_site.py` into `euks.json:literature` and shown under the Eukaryotes table. MAG ids are mapped through the registry (`retired:EEnn (merged into EEmm …)`), so EE14/EE17 show as EE12 and EE29 as EE03; EE13 shows as EE06.

## Outlier-resistant MAG RNA profiles (2026-11-08, user-approved)
- **Why:** a summed MAG profile (reads on all contigs per million) can be set by a few contigs: EE04 94% and EE01 98% of reads on 5 contigs (18S/28S rRNA), EE08 95% (conserved tubulin/ubiquitin/HSP90), EE31 83% (contigs the hull added). That broke the merge rule (EE03+EE31 summed ρ 0.53) and would make spurious network edges.
- **One shared function:** `scripts/robustexpr.py` (`robust_profiles`, `robust_profile`, `summed_profiles`):
  1. per contig and library: reads per kb per million contig-assigned reads (128 Illumina RNA libraries);
  2. keep contigs detected in ≥ 25% of libraries;
  3. scale each kept contig's profile to mean 1;
  4. MAG profile = per-library **10% trimmed mean** of the scaled profiles × the MAG's median contig level.
  - < 5 kept contigs: length-weighted mean (fallback; EE15, EE16). < 20 kept contigs: flagged `thin` (EE08 9, EE10 7, EE12 17, EE36 15): the few dominant contigs can still be a big share of so few.
  - **Deviation from the approved plan (median → 10% trimmed mean):** with sparse eukaryotic contigs the per-library median was 0 in 45% of MAG × library cells (EE02/EE04/EE06/EE07/EE30… had a median profile of 0), leaving no time course and making rank correlations mostly ties at 0. The trimmed mean is 0 in 9% of cells. Outlier resistance comes from the mean-1 scaling (each contig weighs the same) plus the trim. `TRIM = None` restores the median. Merge decisions are the same under either.
- **Used by:** `hull_expr.py` (`out/hull_mags_expr_samples.tsv` is now the robust table; the old summed one is `out/hull_mags_expr_samples_summed.tsv`; `out/hull_mags_expr_method.tsv` has contigs used, method, thin, top-5 read share and ρ robust vs summed), `merge_mags.py` (RNA ρ), `heatmap.py`, `hull_mags.py` (RNA differential abundance), `common.hull_mag_expression()` (Eukaryotes page time courses and medians, Viruses host candidates), and the network node profiles.
  - Status (`network_ok`, ≥ 10 reads per million) and the "expressed neighbour" caution still use the summed `expr_ppm_median`.
  - The split guard (`eukshared.split_parts`) still uses the summed profile; E03's lobes pass either way (it was not rerun: that would need final → … → register_blobs).
- **Robust vs summed (Spearman ρ of the two profiles of one MAG):** EE01 0.01, EE03 0.73, EE04 0.95, EE08 0.76, EE10 0.86, EE12 0.82, EE34 0.72, EE36 0.51; the clean MAGs 0.9–0.99 (EE05 0.98, EE30 0.98, EE23 0.99).
- **Merge rule rerun with robust RNA:** EE03+EE31 ρ 0.996 (summed 0.53), DNA 0.837, AF 2.8/5.5%, BUSCO 55.8 [2.3] + 11.6 [0] → 64.3 [3.9]: **merged** (`retired:EE31 (merged into EE03, 2026-11-08)`); EE29 then re-tested against EE03+EE31: ρ 0.995, DNA 0.859, BUSCO → 68.2 [3.9]. EE06+EE13 still passes (ρ 0.986, DNA 0.827, BUSCO 12.4 + 0.8 → 13.2 [0]). E21 pieces: best ρ 0.87 (EE33/EE35), none pass. EE03 is now E03+E29+E31, 54.0 Mbp, BUSCO 68.2% [D 3.9%], Attheya (genus, medium).
- **RNA differential abundance** now runs on the robust profile. Calls that changed (q < 0.05 gained or lost) are mostly in thin or rRNA-heavy MAGs: e.g. EE12 loses its diesel-vs-control calls and gains ice ÷ water; EE32 gains 2025-26 diesel ÷ control (+2.8 log₂).
- **Cautions:** hulls whose note says `RECHECKED` (E30, E12) are exempt from both automatic hull cautions ("hull added X Mbp" and "hull may include an expressed neighbour"). EE03 now carries "hull added 21.9 Mbp to the 39.4 Mbp source group" (E31 and the E29 hull area are outside the E03+E29 source groups).
- **Data sources on the page:** the time-course charts are titled "Illumina RNA (metatranscriptomes), robust MAG profile, reads per kb per million"; DNA columns say "nanopore DNA"; the heatmap title names both sources.

## Illumina DNA (metagenomes) (2026-11-09)
- **Data:** danaSeq illumina_rna in DNA mode (`--rrna_method none`, human removal, strobealign, alignments covering ≥ 80% of the read) of the CMO Illumina metagenomes against the same coassembly (Grex `~/scratch/project_rnaseq/marine_mg`, local `data_mg/`). Loader: `common.illumina_dna_samples()` / `illumina_dna_contig_matrix()` (lanes summed per library, cached in `data_mg/contig_counts_lib.npz`); `eukshared._mgdna()` / `mgdna_profile()`; `robustexpr.robust_profiles(groups, data='mgdna')`.
- **Libraries:** 169 in the run, 162 used. The 7 dropped are all QC-flagged (`excluded('metagenome')`): low_yield 1343765_plate9_A05, 2599682_plate1_A06, A09, A12 (the four < 0.2 M raw-read failures), pool_mismatch 2599682_plate1_A05 and C08, contamination 2599682_plate1_E02. The ≥ 1 M pairs-on-contigs floor (`MG_MIN_PAIRS`) removes nothing more (smallest used: 5.3 M). 85 of the 162 are pool libraries (CTL/EXP, seawater or sea ice, with a day); the rest are Churchill River incubations, MOS cultures, sites, sea foam and enrichments.
- **Metric (`scripts/mgdna.py`):** length-weighted mean, read pairs on the MAG's contigs per kb per million contig-assigned pairs (reads per kb per million; = mean depth as a share of the library, the nanopore depth × length share divided by length). Not the robust RNA profile: DNA is not spiky like RNA (median top-5-contig share of Illumina DNA reads 9%; EE04 6% against 94% of its RNA); the robust profile on the same DNA agrees (median ρ 0.94) and the length-weighted mean agrees slightly better with nanopore (median ρ 0.58 vs 0.53). MAGs with top-5 share ≥ 0.3 are flagged `few_contigs` (EE01, EE08, EE16, EE33, EE35).
- **Agreement with nanopore DNA** (Spearman per MAG, 36 samples with both, matched by winter + corrected sample name): median 0.58; ≥ 0.8 for EE03 0.94, EE04 0.92, EE07 0.90, EE22 0.87, EE06 0.86, EE27 0.85, EE30 0.84, EE05 0.83; lowest for rare or few-contig MAGs (EE15 0.04, EE36 0.16, EE08 0.17, EE16 0.24, EE35 0.29). P01 0.72.
- **Differential abundance:** `dacore.run_mg` (the same contrasts; ice vs water from same winter/pool/day pairs, 19 pairs), data 'Illumina DNA', `out/hull_mags_da_mgdna.tsv`; `common.hull_mag_da()` concatenates it. Unlike nanopore it has a 2025-26 diesel-vs-control test (12 diesel vs 8 control libraries) and a BACI.
- **Heatmap:** third block, Illumina DNA of the 85 pool libraries.
- **EE30 check (report only, 2026-11-09):** the 98 contigs holding the copies of EE30's 48 duplicated BUSCOs (`out/busco/ee30_dups.tsv`) do not form two coverage clusters in Illumina DNA (2-cluster average linkage on log profiles: 92 / 6). EE30 is thin in the metagenomes (median 0.012 reads per kb per million; a copy contig is detected in a median 42% of pool libraries), so the test is weak. Two things point the way of the "≥ 2 species" call: no duplicated BUSCO has both copies in the old E30 core (30 core contigs, all copy pairs are core–added or added–added), and the copies of one BUSCO correlate less with each other (median r 0.34–0.38) than random pairs of the same contigs (0.39–0.49; Mann–Whitney p 1e-5), while core contigs correlate tightly with each other (0.82). So the second genome most likely sits in the 34.8 Mbp the hull added, with a coverage profile close to but not identical with the core's; Illumina DNA cannot separate it cleanly.
