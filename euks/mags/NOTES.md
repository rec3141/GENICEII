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
3. `hull_panels.py` → **`hulls.tsv`** (curated, hand-picked; never regenerated) → `hull_mags.py` → `plot_hulls.py`.
4. `organelle.py`: exploratory, see Pending. Do NOT run it as part of the pipeline: it rewrites the organelle tables with a flawed split.
- Always run `final.py` before `finalize.py`. `finalize.py` rewrites `contig_euk_group.tsv` from bin groups, and `register_blobs.py` appends the kernel blobs afterwards.

### Rules (see sections above for details)
- **Eukaryote groups (E##):**
  - Built by containment + density gaps (DBSCAN eps 1.5).
  - Seeds must be ≥ 50% Tiara eukarya.
  - GC/24-unit merge; footprint merge.
  - Split rule (`eukshared.split_parts`): Ashman D > 3, > 5 units apart, minor part ≥ 2 Mbp and ≥ 10%.
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
