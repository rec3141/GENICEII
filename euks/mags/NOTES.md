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
