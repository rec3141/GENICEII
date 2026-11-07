# Porticoccaceae 50-400-T64, comebin_2576 (2026-11-07)

This feeds the site section `#porticocc` (site/porticocc.js, site/porticocc.json; suggested nav label "50-400-T64").

## Scripts

- **scripts/extract_cummo.py** writes out/cummo_query.faa.
  - It gets the CuMMO C/A/B proteins for comebin_2576 (contig_0004457) and Cycloclasticus metabat_9342 (contig_0001426) by translating ORFs.
  - The local genemap has no coordinates, so ORFs are matched by exact gene length. The comebin_2576 A protein matches hgt/grex/hc_prot.faa.
- **scripts/porticocc.py** writes out/genes_expr.tsv, out/toolkit.tsv and out/gene_cpm.pkl.
  - Per-gene CPM (of CDS reads), using medians and the within-genome per-mille share.
  - The pooled share (`pshare_*`) is robust when the genome is rare.
  - The toolkit classes are product/KO regexes.
- **scripts/cummo_blast.py** writes out/cummo_blast.tsv and out/cummo_blast.json.
  - Local blastp of subunit A against out/cummo_refs_A.faa: 19 labelled references from NCBI/UniProt, fetched 2026-11-07.
  - An NCBI nr BLAST was submitted but stayed WAITING for more than 50 minutes, so the local set is used. If nr results are saved as out/nr_*.xml, they replace the local table.
- **scripts/porticocc_site.py** writes site/porticocc.json.

## Findings

### Identity
- GTDB: g__50-400-T64 (Porticoccaceae), no species. The genome is 3.4 Mbp, 96% complete and 0% redundant.
- The genus is named after MAG 50_400_T64 (GCA_002162755.1, LBNL), from the Deepwater Horizon plume simulation (Hu et al. 2017 PNAS, doi:10.1073/pnas.1703424114).
- Against that MAG: blastn ANIb ≈ 80.8% over 60% of comebin_2576 (1,020-bp fragments). That makes them the same genus but different species.
- The type MAG also has the CuMMO (A/B/C at 81/79/79% aa) and several cyclohexanone monooxygenases.

### CuMMO (subunit A)
- 81% to the 50_400_T64 type MAG and 75% to another Porticoccaceae MAG (TAL02096).
- 45–48% to the Cycloclasticus short-chain-alkane CuMMOs (Bathymodiolus symbiont; DWH 44_32_T64 and 46_120_T64).
- 38–42% to pmoA/amoA, 24% to archaeal amoA. Actinobacterial references give no significant hit.
- So it is a separate Porticoccaceae CuMMO lineage, probably a short-chain or sub-terminal alkane monooxygenase.
- The Cycloclasticus metabat_9342 CuMMO sits ~48–50% from both the Cycloclasticus alkane CuMMOs and gamma pmoA, so its placement is unresolved.
- The two study CuMMOs share 37% identity with each other.
- A quick FastTree on 21 subunit A sequences was unstable (long-branch artefacts) and was not kept.
- Operon order is C-A-B. Neighbours:
  - upstream: the urea ABC transporter (urease just beyond), CoxE/MoxR and an acyl-CoA dehydrogenase;
  - downstream: a multicopper oxidase, a ResB-like protein and an RND transporter.

### Toolkit
- Absent: no alkB, no CYP153-type annotation and no ring-hydroxylating dioxygenase.
- Present:
  - 3 cytochrome P450s;
  - 13 Baeyer–Villiger MOs (cyclohexanone/acetone MO, K03379/K18371/K18091/K21730);
  - cyclopentanol dehydrogenase (K17754);
  - 12 alcohol and 9 aldehyde dehydrogenases;
  - 86 β-oxidation genes;
  - methylcitrate cycle;
  - glyoxylate shunt.
- Isocitrate lyase is its top transcript, about 29‰ of the genome's transcripts (median, 2025-26 diesel).
- Pooled within-genome shares (2025-26 diesel vs background = control plus both pools pre-dosing) are flat for every group: |log2| ≤ 0.5. CuMMO is 5.0 vs 6.9‰.
- So the oxidation programme is constitutive. Diesel selects for the organism rather than inducing the genes.

### Ecology
- 2025-26, before dosing: present in both pools (EXP 300–1,500 TPM, CTL ~250).
- It rises after Cycloclasticus lorbin_184:
  - lorbin_184 is ≈20k TPM at days 6–8 and plateaus at 120–165k from day 17;
  - comebin_2576 is 4k at day 6, 27k at day 20 and 70k (7.1%) at day 82, still rising.
- Diesel vs control after dosing: 35,915 vs 282 TPM. BACI +3.6, q 0.02.
- 2024-25: ~7k TPM at days 55–60, behind metabat_9342. BACI +3.9, q 0.68.
- DNA:
  - Illumina k-mer share is 0.002–0.02% before dosing and 0.11% by day 13 (2025-26). In 2024-25 it is 0.2–0.35% at days 55–60, with the control below 0.02%.
  - Nanopore share is 0.01–0.06% before dosing and 0.06–0.17% at days 4–6, reaching 0.19% in the DAY04 oil-snow sample.
  - The June 2025 surface sample of the 2024-25 diesel pool is ~4× the control's.
- Units: out/dna_rel_abundance.tsv is already in percent. The "~8% of nanopore DNA" figure in the brief was a ×100 slip; the real value is ~0.08–0.2%.
- Ice vs water: CLR −0.89, lower in ice in 20/24 pairs (q 0.03). Top vs bottom ice: −0.87 (4/13 pairs up, not significant).

### Mobile elements
- The CuMMO contig is chromosomal: no island, plasmid or provirus.
- 2 IslandPath islands (GI_1486 14 kb, GI_2506 8 kb), with no hydrocarbon genes.
- 3 geNomad plasmid-like contigs, 1 of them an ICE (contig_0064327, with tra/virB4).
