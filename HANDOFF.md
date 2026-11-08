# Handoff: GENICE II CMO diesel exploration (state at 2026-11-08)

Read this first when resuming. The repo README has the per-analysis sections; each analysis folder has a NOTES.md.

## Where things live
- **Analysis repo:** `/data/dev/rnaseq-exploration`, git, `main` → https://github.com/rec3141/GENICEII (public).
  - `.gitignore` is a whitelist: code, md and `site/**` are tracked; small tables go in with `git add -f` (< 1 MB).
- **Live site:** https://cryomics.org/GENICEII/
  - Deploy only with `scripts/deploy_site.sh` (deploys `origin/main` and stamps `?v=<sha>` on the script tags).
  - Never deploy uncommitted work.
- **Python:** `/data/miniforge3/bin/python3`.
- **Headless render checks:** `google-chrome --headless=new … --dump-dom/--screenshot` against `python -m http.server <port>` started with `timeout`. Never `pkill -f` (it kills its own shell).
- **Metadata archive:** `/data/dev/genice-sample-metadata` (not a git repo).
  - QC flags live in `sources/sample_corrections.tsv`; backup `.bak-20261108`.
  - `scripts/common.py` `EXCLUDE_FLAGS = (low_yield, pool_mismatch, contamination)` drops flagged libraries everywhere: 128 RNA libraries are used.
- **Grex:** reach it via `ssh -o BatchMode=yes concentration 'ssh -o BatchMode=yes grex "…"'`.
  - One call at a time. No loops without `ssh -O check`. Never retry on exit 255.
  - `mag_marine` is read-only.

## Data state
- **Site runs on FULL-DEPTH RNA** (Grex `~/scratch/project_rnaseq/marine_full`, 4.37 G pairs). The capped inputs are kept in `data_capped/`.
- **Gene table:** the old Bakta tags with coordinate-corrected KOs. Feature order is identical to the full-depth run.
- **Excluded RNA libraries:** both day-18 control libraries (`V350357237_L01_D09`, `L02_D10`) and their metagenome (`2599682_plate1_E02`), flagged `contamination` for a transient Zhongshania/Oleispira hydrocarbon-degrader community. Merged in PR #16.
- **Bin domains** are rule-based. `scripts/bin_domains.py` writes `data/bin_domain.tsv` (Tiara length-weighted: ≥ 50% eukarya → eukaryote; ≥ 50% organelle → organelle).
  - Helpers in `common.py`: `special_label`, `is_prokaryote`, `is_plastid`.
  - The user wants rules, not one-offs, for reproducibility.
- **Lineage subtitles:** `scripts/lineage.py` writes `site/lineage.json`; `fig()` in `site/app.js` adds a subtitle when a chart title names exactly one MAG (`MAG_RX`).

## Merged PRs
#1–#16, all deployed:
- the pages (nitrogen, 50-400-T64, oxygen, eukaryotes, states, light + under-ice model, guilds);
- dropdown nav and theme picker;
- bin domains and lineage subtitles;
- the full-depth rebuild;
- the day-18 QC exclusion.

## OPEN / IN PROGRESS
1. **PR #17** (branch `euk-groups`, local checkout is on this branch with uncommitted changes in `data/` and `euks/mags/`). It wires the eukaryotic genome groups into the domain labels and pages:
   - `data/euk_groups.tsv`, `data/contig_euk_group.tsv.gz`;
   - the `euk_group` column in `bin_domain.tsv`;
   - in `common.py`: `euk_group`, `euk_group_info`, `is_euk_group`, `euk_group_label`, `contig_host_label`;
   - the Eukaryotes page rewrite, the Viruses page (eukaryotic host candidates; E07 hosts giant virus `contig_0000564` ρ 0.87 and virophage `contig_0049668` ρ 0.74), and E## in `lineage.json`.

   A background task was refreshing it with:
   - the E03 split (E03 24.8 Mbp + new E31 7.6 Mbp, rule `split_parts` in `euks/mags/scripts/eukshared.py`);
   - kernel-blob groups E25–E30;
   - virus groups V01 (= old E09) and V02 (`contig_virus_group.tsv`, label "Virus V02");
   - DA columns from `da_blobs_wide.tsv`;
   - `MAG_RX` extended to `V\d\d`.

   **Do not merge or deploy without the user's go-ahead.** The user usually says "merge and deploy".
2. **EE##/VV## hull-curated MAGs** (the user's "final refinement"). Eyes draw hulls; data records them:
   - Per-group zoomed t-SNE panels with a grid go in `euks/mags/out/hulls/panels/`.
   - Polygon vertices are picked by eye from the panels and recorded in **`euks/mags/hulls.tsv`** (id, source_group, vertex_order, x, y, note). This is curated input, never regenerated.
   - `euks/mags/scripts/hull_mags.py` does point-in-polygon assignment of any contig:
     - eukaryotic hulls exclude Tiara bacteria/archaea;
     - overlaps go to the source group, otherwise the nearest centroid.
   - IDs: EE## = source E number (E05 → EE05), VV## from V##; recorded in `out/id_registry.tsv`.
   - Outputs: `out/hull_mags.tsv`, `out/contig_hull_mag.tsv`, `out/euk_hull_mags_tsne.png` (3 panels with hull outlines) plus verification panels.
   - Then wire EE/VV into PR #17 the same way as E##/V## (data/ copies, `common.py` helpers, Eukaryotes/Viruses pages, `lineage.json`).

   Check `euks/mags/NOTES.md` and whether `hulls.tsv` exists to see how far it got.

   **Done 2026-11-08 (PR #17):** EE## and P01 are wired as the site's labels (Eukaryotes table, time courses, figure `site/img/euk_hull_mags_tsne.png`; host labels on Light, Nitrogen, Sulfur; Viruses host links; `lineage.json`; `MAG_RX` matches EE/VV/P/E/V). E##/V## are provenance only.
   - VV01 was retired: it is a Caudoviricetes cluster, not a genome. Viral genome-group rule: ≥ 50% of bp in giant/eukaryotic-virus classes. There are no VV## now.
   - The Viruses page lost its "Large viral genome groups" section; one line covers the phage cluster and the former V02 (now P01).
   - Details: `euks/mags/NOTES.md`, "Site wiring of the hull MAGs".
3. **Food-web network (next, requested):** MAG–MAG network analyses on the Food web page.
   - Nodes: prokaryotic MAGs plus network-ready eukaryotic MAGs (`data/hull_mags.tsv` `network_ok`, 21 EE; mind the `caution` column, e.g. EE01). There are no viral genome groups (VV01 retired), so viruses enter as vOTUs.
   - Ideas already discussed: co-response guilds (`guilds/` page) and the genome pairs flagged there.
3b. **Organelle groups (user request 2026-11-08):** groups ≥50% Tiara organelle become P## (plastid) or M## (mito) and are excluded from the E/V series. The old V02 was a plastid that geNomad mis-called viral. They are saved for later; there are no hulls and no site wiring yet.
   - Pending: 486 plastid genes means several plastid genomes merged together. Next steps:
     - count single-copy markers (psbA/rbcL/psaA…) per P group;
     - split by coverage (`split_parts` rule);
     - match each plastid to a host E## by ρ and check the lineages agree;
     - check that all coassembly organelle contigs fall inside the P/M groups.
3c. **Viruses:** VV01 was a cluster of tailed-phage contigs (2,165 contigs, 21 Mbp, Caudoviricetes) and has been retired. V/VV genome groups now require giant-virus (Nucleocytoviricota) content. A dedicated danaSeq virus module is planned (user, 2026-11-08), so don't build deeper virus binning here. The jumbo phage `contig_0000248` (569 kb, 405×) is worth a look later.
4. **Optional / not requested yet:**
   - IBP host labels in `ice/cold/scripts/ibp.py` still use Kaiju, not the eukaryote groups.
   - Put the metadata archive under git.

## euks/mags pipeline
- **Order:** survey → compact → final → expr → finalize → kernel → da → register_blobs → plot, then hull_mags once `hulls.tsv` exists.
- **Candidate bins:** the config `LIST` at the top of `survey.py`.
- **Statuses:** network / "little CMO RNA" (pure but mostly non-CMO reads; not discounted as genomes) / doubtful.
- **danaSeq issue #129** (https://github.com/rec3141/danaSeq/issues/129) covers why eukaryotes are lost in `mag_analysis`: the Binette 10 Mbp cap, DAS Tool, no eukaryotic QC/taxonomy.

## Full-depth rerun on Grex (done)
- **Driver:** `~/scratch/project_rnaseq/rna-marine-full-grex.sh` (OUT=`marine_full`, `MAX_PAIRS=0`).
- **Jobs:** 8159221–28, all COMPLETED.

## Conventions the user expects
- **Commit and PR trailers:** commits end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`; PR bodies end with `🤖 Generated with [Claude Code](https://claude.com/claude-code)`.
- **Git workflow:** work on a branch → open a PR → merge only on the user's say-so → `./scripts/deploy_site.sh`.
- **No one-off installs** on Grex or locally; use existing envs or containers.
- **Plain-language page text.** Season is a bundle of drivers (light, temperature, ice, nutrients), not a nuisance variable.
