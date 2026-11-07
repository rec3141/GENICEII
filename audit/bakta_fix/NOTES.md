# Check of the BAKTA_EXTRA locus-tag fix (danaSeq f2fc406, PR #100), 2026-11-06

- On Grex, mag_marine/annotation/bakta/extra was rewritten in place on 2026-11-02 (locus_tag_map.tsv present). The marine illumina_rna summary was rerun on 2026-11-04.
- This analysis (site up to v11) was built 2026-10-21..26 on the OLD extra-run tags. Every join went either through the coordinate KO remap (data/ko_remap.tsv.gz) or through files taken from the extra run itself (CANT-HYD, Pfam cold search, amino-acid counts, island features). The site is therefore unaffected.
- One stray use of the pipeline's misaligned ko_tpm, in fig/f4_hc_genes.png (K00217 line), was switched to data/ko_cpm.tsv.
- Pipeline ko_tpm (Nov 4, ko_tpm_20261104.tsv here) vs the local corrected KO CPM, log1p mean over 130 samples:
  - new r = 0.965; old r = 0.848;
  - K00217: old 17.2 → new 2.5 TPM (local 3.3 CPM);
  - K16242 (dmpN): 2.8 → 37.6 (local 46.2);
  - K00496 (alkB): 127 → 154 (local 226).
- The remaining differences are TPM (length-normalised) vs CPM of read pairs.
- Local gene IDs (genemap.pkl, gene_counts.npz) still use the old tags. To join with new Grex outputs, translate through mag_marine/annotation/bakta/extra/locus_tag_map.tsv.
