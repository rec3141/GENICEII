# Eukaryotes (2026-11-07)

Site section `#euks` (site/euks.js, site/euks.json). Nav label: "Eukaryotes". Code:

- **scripts/euks.py.** Main analysis.
  - Composition: 18S from SILVA rRNA read taxonomy (genus rank partitions each domain).
    - 21 groups by lineage regex, each with a trophic role.
    - Likely contaminants removed (mammals, fish, insects, land plants).
    - Shares are % of eukaryotic 18S.
  - Activity:
    - Kaiju mRNA has no eukaryote database, so plastid activity is measured as Kaiju "Cyanobacteriota" mRNA % (plastid-encoded transcripts) plus plastid MAG TPM.
    - Viruses: Megaviricetes and Maveriviricetes CPM from viruses/out/vrna_samples.tsv; virophage contig_0049668 from votu_rna_cpm.pkl.
  - Contrasts:
    - Seawater diesel vs control after dosing, per winter.
    - BACI per winter.
    - Ice vs water (24 pairs) and top vs bottom (13 pairs), via ice/scripts/contrasts.py.
  - Correlations: Spearman of each genus/group with giant-virus and virophage CPM, all pool samples per winter, BH-adjusted.
- **scripts/euks_site.py.** Builds site/euks.json.
  - Includes a plastid-bin table and grazer-vs-Cycloclasticus correlations (2025-26 diesel pool, SW+SI).

## Findings (seawater medians after dosing; one pool per treatment, descriptive)

**Algae collapse in the 2025-26 diesel pool.**
- Plastid-encoded mRNA: 1.8% vs 17.6% (control); pre-dosing 17.8% vs 19.1%; BACI −3.1 log2.
- Plastid MAG TPM: BACI −2.4.
- Group shares of eukaryotic 18S, diesel vs control (BACI in brackets):

| Group | Diesel | Control | BACI |
|---|---|---|---|
| Haptophytes | 0.04% | 3.0% | −3.5 |
| Dinoflagellates | 0.17% | 1.7% | −2.3 |
| Cryptophytes | 1.1% | 7.6% | −1.9 |
| Diatoms | 11% | 36% | −1.7 |

- Thalassiosira, Fragilariopsis, Micromonas and Bathycoccus are all down.

**Bacterivorous and mixotrophic flagellates replace them (2025-26).**
- Chrysophytes and synurophytes reach 58% vs 11% (BACI +2.2): Synura, Mallomonas, Paraphysomonas, Chrysolepidomonas.
- Choanoflagellates: 5.8% vs 1.4% (BACI +1.5).
- Within diesel-pool samples, these track Cycloclasticus TPM: telonemids ρ 0.69, bicosoecids 0.58, MAST 0.44, Cercozoa 0.43, chrysophytes 0.41.
- These fall as Cycloclasticus rises: diatoms −0.60, haptophytes −0.61, dinoflagellates −0.79, katablepharids −0.80.
- Time since dosing is a confounder.

**2024-25.**
- Plastid mRNA: 4.9% vs 14.8% (pre-dosing 14.5% vs 12.3%).
- Higher in diesel: Cercozoa +0.9 and katablepharids +3.1 log2.
- Chrysophytes did not bloom.
- Weak baseline: single pre-dosing sample per pool, and the pools were filled two weeks apart.

**Giant viruses.**
- 2025-26 Megaviricetes: 4,639 CPM pre-dosing → 64 in the diesel pool, while the control is at 2,358 (BACI −4.5).
- Correlated with Mamiellophyceae (Micromonas ρ 0.69, Bathycoccus 0.65 in 2024-25; Bathycoccus 0.78 in 2025-26), Thalassiosira (0.82) and katablepharids (0.85).
- Interpretation: they follow their hosts (prasinoviruses etc.).

**Virophage contig_0049668 (2025-26, control pool only).**
- Correlates best with Phaeocystis antarctica (ρ 0.80, q 1e-10), haptophytes as a group (0.79) and Pedinellales (0.79).
- Bicosoecids/Cafeteria are negatively correlated (−0.45).
- Hypothesis: a Phaeocystis–giant virus–virophage system, absent from the diesel pool because haptophytes were wiped out.

**Ice vs water.**
- Higher in ice:
  - Fungi: 23 of 24 pairs (+1.6; Helotiales, Capnodiales).
  - Cryothecomonas: 23 of 24 (sea-ice diatom parasitoid).
  - Protaspis: 23 of 24.
- Lower in ice: kinetoplastids (higher in only 3 of 24 pairs).

**Top vs bottom ice.**
- Green algae are higher in top ice in 12 of 13 pairs.
- Dinoflagellates are higher in top ice in 0 of 13 pairs.

**Caveats.**
- 18S share mixes cells with copy number.
- It shows composition, not biomass.
- Plastid mRNA % is diluted by the bacterial bloom (≤ ~2×).

## Bin domains (2026-11-07)
Labels come from `special_label`, using the Tiara-based `data/bin_domain.tsv` (see the README). Eukaryotic and organelle bins are relabelled, and they are excluded wherever a bacterial genome is assumed. The three eukaryotic bins are comebin_4527, comebin_3937 and lorbin_494_sub.
