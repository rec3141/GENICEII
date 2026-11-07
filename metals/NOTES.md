# Metals & toxicity (2026-11-06)

No metals or trace elements were measured: the CMO 2024/2025 logs only note that nutrient samples were taken. Diesel is low in metals (unlike crude: V, Ni porphyrins). Everything here is inferred from gene content and transcription.

## Method
- **Gene sets** (`scripts/msets.py` → `out/members.tsv`):
  - 18 sets from the corrected KOs (genemap.pkl) plus Bakta product regexes. A gene joins the first set it matches.
  - Sets cover toxic-metal resistance (Hg, As, Cu, heavy-metal RND, Zn/Cd/Pb, chromate, tellurite, Ni/Co export, metallothionein), micronutrient homeostasis (Fe storage/Fur, Fe uptake, Mn/Zn uptake, Ni/Co uptake), solvent stress (RND multidrug efflux, cti, psp), oxidative stress and protein-folding stress.
  - The confident tier drops KO-only hits whose product is a bare hypothetical/DUF.
- **Expression** (`scripts/mexpr.py`):
  - Units: CPM of CDS read pairs, with corrected labels.
  - Diesel vs control: medians after dosing, plus BACI against each pool's pre-dosing seawater. In 2024-25 the diesel pool has one pre-dosing library after the day-5 correction.
  - Ice vs water and top vs bottom: the paired contrasts from ice/scripts/contrasts.py.
  - Hosts: category × bin CPM.
  - Within-host: category reads per million host CDS reads, over libraries with ≥5,000 host reads.
- **Mobility** (`scripts/mmob.py`):
  - Same flags as hgt/scripts/s3_context.py, plus co-location with confident hydrocarbon genes within 10 kb.
  - Coordinates come from the reconciled extra Bakta run on Grex (`grex/metal_ctx.sh`, job 8157239), translated back to the pre-Nov-2 tags via locus_tag_map. All 143,645 members got coordinates.
- **Site:** `scripts/mbuild.py` → site/metals.json; page in site/metals.js.

## Findings
- **Diesel, 2025-26 seawater, BACI (log2):**
  - psp +2.97, Zn/Cd/Pb export +2.35, cti +2.22 (tiny base: 2 CPM), metallothionein +1.62 (1 CPM), arsenic +1.33, Fe storage +1.25, copper +1.05, heavy-metal RND +0.99, mercury +1.11.
  - Raw diesel/control: psp +2.36, cti +2.22, Zn/Cd/Pb +1.38, arsenic +1.38, Fe storage +1.22.
- **Diesel, 2024-25 BACI** (one pre-dosing library): copper +2.17, cti +1.85, Ni/Co uptake +1.70, Ni/Co export +1.26, mercury +1.24, Zn/Cd/Pb +1.21. psp −0.63.
- **The community-level metal rise is mostly the Cycloclasticus bloom.**
  - In 2025-26 diesel seawater, Cycloclasticus lorbin_184 carries copper 521 of 3,140 CPM (control 11), Zn/Cd/Pb 75 of 335, arsenic 59 of 489 and mercury 46 of 235.
  - Its within-host levels are not higher after dosing: copper 4,148 vs 6,097 before dosing; Zn/Cd/Pb 852 vs 805 per million host transcripts.
  - Sulfurimonas semibin_088 carries arsenic 114 CPM (diesel) vs 8 (control).
- **psp is a real within-host response of UBA9145 vamb_4149 (Pseudohongiellaceae).**
  - 2025-26: 43,047 per million host transcripts in the diesel pool vs 2,177 before dosing (n = 25 vs 5) and 2,557 in the control.
  - 2024-25: 19,400 vs 6,030 before dosing (n = 1) and 1,473 in the control.
  - It contributes 1,038 of 1,406 CPM of the community psp in 2025-26 diesel seawater.
- **Oxidative stress within hosts:**
  - Oleispira lorbin_369: 6,634 in the 2024-25 diesel pool vs 721 before dosing (n = 1) and 737 in the control; 15,918 in the 2025-26 diesel pool.
  - Cognaticolwellia lorbin_135: 28,069 in the 2025-26 diesel pool vs 7,942 in the control (n = 1).
  - UBA9145: 30,336 in the 2024-25 diesel pool vs 7,111 before dosing.
- **Ice (2024-25, top − bottom; 13 pairs):**
  - Higher at the top in 13/13 pairs: mercury (+0.55), arsenic (+0.42), Fe storage (+0.54).
  - 12/13: copper (+0.57), tellurite (+0.41).
  - Heavy-metal RND 11/13, oxidative stress 10/13.
  - Ice vs water (24 pairs) is weak for metals: |delta| < 0.3 and q > 0.1. psp is +0.57, 17/24.
  - This is consistent with brine concentration at the top (61–166 PSU), but composition is a confounder.
- **Mobility** (share mobile-associated; strict-plasmid share in brackets):
  - Ni/Co export 20.6% (2.5%), mercury 18.5% (3.8%), Zn/Cd/Pb 15.7% (4.4%), copper 12.8% (3.4%), heavy-metal RND 12.3% (3.0%).
  - Hydrocarbon genes (confident): 6.8% (1.3%).
  - psp is 17% mobile, mostly islands (9.1%).
- **Co-selection:**
  - 12 of the 51 strict geNomad plasmids with confident hydrocarbon genes also carry toxic-metal genes; 214 strict plasmids carry metal genes in all.
  - Examples: Boseongicola lorbin_194 (442 kb; As, Cu, Zn/Cd/Pb + RHO dioxygenase, meta-cleavage), unclassified vamb_488_sub (421 kb), Sulfitobacter metabat_1122 (337 kb; As, czc + intradiol), Hoeflea comebin_1575.
  - The Sphingorhabdus CYP153/alkB plasmid (contig_0009679) has no metal genes. The Pseudosulfitobacter alkB plasmid-like contig_0009104 has one czcA.
  - Only about 1–2% of metal genes sit within 10 kb of a hydrocarbon gene.

## Caveats
- One pool per treatment, and the pools hold different strains (Ice vs water popgen).
- 2024-25 BACI rests on one pre-dosing library.
- Broad families (TonB receptors, RND, MerR) overlap non-metal functions.
- KO-only members on hypotheticals are excluded from the headline numbers but kept in out/cat_cpm_all.tsv.
