# Cycloclasticus genomics in the CMO diesel mesocosms (2026-10-23)

The EXP pool got diesel in both seasons: 2024-25 on DAY05 at 16:30, 2025-26 about 136 L on DAY01. CTL is the unoiled pool.
Data comes from the marine Flye coassembly in mag_marine (DAS Tool bins, Bakta), the 130 BGI metatranscriptomes and 72 CMO nanopore metagenomes.
Website data is in `cyclo.json`. Scripts are in `scripts/`.

## Genomes
Five DAS Tool bins are Cycloclasticus, falling into three species-level groups:

| MAG | size | contigs | N50 | GC | CheckM2 compl/contam | closest genome (skani, GTDB r232 / NCBI) |
|---|---|---|---|---|---|---|
| metabat_9342 | 2.78 Mb | 21 | 265 kb | 44.5% | 99.99 / 0.8 | **98.7% ANI** to sp047634535, a MAG from a Central Arctic Ocean sea-ice core (87°N, MOSAiC 2020); 94.6% to the Kivalliq (Nunavut) sediment + crude-oil MAG |
| lorbin_184 | 2.67 Mb | 46 | 92 kb | 42.0% | 94.7 / 10.2 | 97.7% to sp001566825 (Canadian marine MAG, 2013) |
| metabat_5432 | 2.87 Mb | 65 | 190 kb | 42.6% | 100 / 12.8 | 96.7% to sp019090295 (Douglas Channel BC, dilbit-coated beads, 360 m) |
| comebin_3787_sub | 1.62 Mb | 118 | 20 kb | 41.6% | SCG 67% | 95.8% to lorbin_184, 93.7% to sp001566825 (sibling strain or species) |
| lorbin_722_sub | 0.64 Mb | 52 | 15 kb | 41.8% | SCG 53% | 98.0% to sp001566825, 94.7% to lorbin_184 (fragment of a sibling strain) |

- None of the MAGs is within skani range (<~80% ANI) of the cultured isolates. C. pugetii PS-1/PY97M, C. zancles 78-ME, sp. P1 and PY97N form one species cluster at 98.5–99% ANI with each other.
- AAI to C. zancles is 82% (lorbin_184) and 74% (metabat_9342). These are uncultured, cold-water and sea-ice lineages, whose closest relatives all come from Canadian or Arctic oil or ice settings.
- No contig is circular. lorbin_184 and metabat_5432 carry about 10–13% CheckM2 contamination.
- More Cycloclasticus sits outside the final bin set:
  - 1,162 unbinned contigs are classed Cycloclasticus by Kaiju, carrying 123 hydrocarbon genes.
  - Several non-selected bins match C. pugetii (comebin_3882) and sp033791505 (comebin_3580, semibin_637, vamb_076).

## Hydrocarbon gene inventory (Bakta products; CANT-HYD at noise cutoffs)
- **PAH and aromatic ring-hydroxylating dioxygenases.** Each genome has 14–21 alpha subunits and 12–20 beta subunits. CANT-HYD counts NdoB 6–8 and NdoC 2–5 per complete MAG, and MAH alpha/beta 1–3. This multi-copy RHO repertoire is the Cycloclasticus signature.
- **Complete aromatic pathways.**
  - Extradiol cleavage: biphenyl-2,3-diol / 1,2-dihydroxynaphthalene dioxygenases and catechol 2,3-dioxygenase.
  - The full meta-cleavage lower pathway, through 2-hydroxymuconate semialdehyde dehydrogenase/hydrolase, 2-oxopent-4-enoate hydratase, 4-hydroxy-2-oxovalerate aldolase and acetaldehyde dehydrogenase.
  - Salicylate hydroxylase, and a phenol hydroxylase (dmpKLMNOP) plus toluene-4-monooxygenase (tmoABCDE) cassette.
- **Alkane genes** (unusual for a PAH specialist):
  - metabat_9342 has **CYP153** (CANT-HYD score 670) and a **Cu-membrane monooxygenase pmoCAB-like operon** (contig_0001426; K10946/K28504/K10945).
  - alkB (K00496, Bakta "Alkane 1-monooxygenase") is in every MAG, but none passes the CANT-HYD AlkB noise cutoff, so these are divergent alkB/FA-hydroxylase-like proteins.
  - The CYP153 and pmo-like genes are shared only with the Kivalliq crude-oil MAG (95–99% aa) and the Arctic sea-ice MAG (CYP153 100%). The isolates lack them.
- **Chemotaxis.** 9–11 chemotaxis/MCP genes in metabat_9342 and metabat_5432.

## Expression
**Who responds, and when**
- **2024-25.** Only metabat_9342 (the sea-ice species) responds at all; the whole 2024-25 bloom is this one species. Its aromatic operons come on at day 7, about 1.5 days after diesel, and rise more than 100× to 1,600–2,200 CPM by days 55–63.
  - The most expressed are the pmoCAB-like operon (pmoC 264 CPM, about 60× control), CYP153 (112 CPM, 80×), and the benzoate-DO/dmp/tmo/meta-cleavage operon (op14/op15, contig_0004615).
  - lorbin_184 and metabat_5432 are absent at DNA and RNA level until late. Cycloclasticus reaches 42–51% of all mapped transcripts at days 55–63, including the marine oil snow samples.
- **2025-26.** The lorbin_184 complex dominates from day 2. Its operons start from a non-zero baseline before oiling (pre 16–74 CPM), roughly double by days 2–10, and reach 250–900 CPM later.
  - Its alkB (JPFHOD_002031207) reaches 108 CPM vs 1.3 in control.
  - metabat_9342 follows from day 6, and metabat_5432 stays minor.

**Induction is mostly population growth, not per-cell regulation**
- Within-genome TPM of lorbin_184 hydrocarbon genes (2025-26, the only season with pre-oil Cycloclasticus RNA) changes little: median log2 FC is −0.75 early vs pre.
- Aromatic genes are constitutively a large share of the Cycloclasticus transcriptome. The community-level increase comes from Cycloclasticus expanding and being transcriptionally active.

**RNA:DNA ratio (nanopore DNA, 2025-26)**
- lorbin_184 sits at about 1 (0.4–1.3) before oil, and rises to 2–4 by days 2–6 in EXP.
- metabat_9342 shows 2–12× by day 6. In 2024-25 it shows 8–85× at days 7–8, where its DNA is still only 0.01–0.03% of reads.
- So the founder species was transcriptionally hyperactive before it bloomed.

## DNA abundance: different species each season
- **2024-25.** All Cycloclasticus are essentially absent in DNA through day 10 (≤0.03% of mapped bases; metabat_9342 is the only trace). In RNA, the 2024-25 bloom is almost entirely **metabat_9342**: 19–51% of all mapped transcripts at days 55–63. lorbin_184 stays ≤0.09% for the whole season.
- **2025-26, before oil.** A different species, **lorbin_184** (plus its siblings comebin_3787_sub and lorbin_722_sub), is already present in **both** pools:
  - lorbin_184: 0.42–0.72% of mapped bases in EXP vs 0.07–0.12% in CTL. Except EXP PRE03, a low-depth outlier at 0.20%.
  - lorbin_184: 0.3–0.7% of transcripts before oil, rising to 10–12% by days 20–82.
  - metabat_5432 shows the reverse, 0.14–0.18% in CTL vs 0.05% in EXP.
  - metabat_9342 is low in 2025-26 (0.01–0.04% before oil), rising to 4.5% by day 82.
- **Not from the intake.** Churchill River site samples from Nov 2025 (SITE1/3/4) have no Cycloclasticus.
- **Interpretation.** The 2025-26 pool populations are resident in the facility, not brought in with the intake. The EXP pool held about 5× more lorbin_184 than CTL before oil. But the species that bloomed in 2024-25 (metabat_9342) is not the one waiting in 2025-26, so this is not a simple carry-over of last year's bloom. Spring/summer 2025 colonisation of the pools, or legacy hydrocarbons selecting a different species, are both possible; there are no samples between Jan and Nov 2025 to tell. The pre-existing population explains the faster 2025-26 response (day 2 vs day 7, and a large bloom by day 11–20 vs day 55).

## Mobile elements and HGT
- **The dmp/tmo/meta-cleavage cassette is not plasmid-borne.**
  - geNomad scores lorbin_184 contig_0030642 (30.7 kb) as a plasmid (0.94), but its ends carry housekeeping genes (fabA/fabB, malic enzyme, flgE), with no mobility genes.
  - It is colinear with C. zancles 78-ME chromosome 2.117–2.142 Mb. Its depth matches the chromosome (copy ratio 1.1, n=9 samples).
  - Read it as a chromosomal aromatic-catabolism island.
- **Evidence for exchange between species.**
  - Cassette proteins are 86–98% identical to C. zancles (mean about 91%), above the 82% genome-wide AAI. dmp/tmo homologs are absent from the other four isolates.
  - metabat_9342 has a rearranged version: tmo + dmp + meta + benzoate-DO on contig_0004615, about 77% nt identity, chromosomal.
  - comebin_3787_sub has the lorbin_184 version at 92–95% nt identity, the same as its genome ANI. That pair is co-inheritance, not recent transfer.
- **Real plasmid, no hydrocarbon genes.** lorbin_184 contig_0003999 (119 kb) is a conjugative plasmid (MOBH relaxase, T4SS virB4/t4cp1, 131 genes) at about 1 copy per chromosome.
- **Other plasmid-like contigs:**
  - Small plasmid-like contigs carrying hydrocarbon genes: metabat_5432 contig_0380606 (6.8 kb, benzene DO alpha/beta, transposase within 5 genes), lorbin_184 contig_0047870/0011601 (plasmid 0.71) and metabat_9342 contig_0967316.
  - 32 hydrocarbon genes on unbinned Cycloclasticus contigs are plasmid-like (geNomad ≥0.7). One example is a tmoABCDE cluster on contig_0184737, plasmid score 0.94.
- **Islands and integrons.** IslandPath islands overlap only 2 hydrocarbon genes (metabat_5432: xylene MO subunit and a dihydrodiol dehydrogenase, both next to mobility genes). Integron attC/elements on Cycloclasticus contigs never contain hydrocarbon genes.
- **Defense.** Only type III-B CRISPR (Cmr5) genes were found. DefenseFinder output in mag_marine is empty.

## Caveats
- The KO joins in the original genemap were misaligned (eggNOG/kofam were run on the Sep-24 Bakta IDs). Fixed by the coordinator's coordinate remap; this analysis uses Bakta products plus the corrected KOs.
- One pool per treatment, so time points are pseudoreplicates. Nanopore depth is low: Cycloclasticus depth is ≤6×, and DNA is 2025-26 days ≤6 plus 2024-25 days ≤11 only.
- CPM is relative to all CDS pairs in the sample. Within-genome TPM is only computed where a MAG has ≥3,000 CDS pairs.

## CANT-HYD at trusted cutoffs (all proteins, job 8119748)
- Across the whole assembly, 490 genes pass --cut_tc. Cycloclasticus TC hits are in `cyclo.json` under `canthyd_tc`; NC hits under `canthyd`.
