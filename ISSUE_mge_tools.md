# mag_analysis MGE tools: DefenseFinder empty, MacSyFinder without loci (mag_marine)

Found on 2026-10-26 while building the genomic-island analysis. No code has been changed yet.

## DefenseFinder: empty output that is reported as success
- **Symptom.** `mag_marine/mge/defensefinder/{systems,genes,hmmer}.tsv` contain only header lines.
  - The task (Sep 24, workDir /data/work/4c/f31b13…) finished in 90 s on about 13 M proteins, with exit 0.
- **The tool itself works.** On Grex (job 8124280, `~/scratch/project_rnaseq/explore/dftest/`), the same container and `bin/parallel_defensefinder.py` were run on one MAG (Marinomonas lorbin_161, 4,326 proteins) with `--models-dir` pointing at the mounted `defensefinder-models.sqsh`. It found 11 systems.
- **Likely cause.** Every defense-finder chunk failed in the original run, and the wrapper swallowed the failures:
  - The wrapper logs a non-zero exit from a chunk only as `[WARNING] Chunk i defense-finder exited …`.
  - It then merges whatever chunk directories exist and writes header-only TSVs.
  - The process script runs it under `set +e` and exits 0.
  - `defensefinder_models` defaults to null, so the run depends on `defense-finder update` succeeding at run time. If the models were missing or unreadable, each chunk would fail this way.
  - The work directory no longer exists, so the exact error can't be recovered.
- **Fix options.**
  1. In `parallel_defensefinder.py`, exit non-zero when every chunk fails, or when more than half of them fail.
     - Alternatively, write a `FAILED` marker that the process turns into a task error.
  2. Require `--defensefinder_models` (the sqsh image) on HPC runs, rather than `defense-finder update`.
  3. Rerun DefenseFinder for mag_marine. To keep mag_marine read-only, run it as a separate job writing to a side directory.

## MacSyFinder (CONJScan/TXSScan): systems without loci
- **What was run.** `--db-type unordered` on the whole-assembly `annotation.faa`, following the comment that unordered mode "is appropriate for fragmented contigs".
- **Why it doesn't help here.**
  - In unordered mode MacSyFinder treats the input as one bag of genes. It reports whether a system's components occur somewhere in the bag, not where.
  - For a 1.55 M-contig coassembly, the 131k reported "systems" (each CONJScan model repeated about 4–6k times) can't be placed on contigs or MAGs.
  - So they can't mark conjugative elements or islands.
- **Fix.** These are long-read contigs, so use `--db-type gembase` with per-contig replicons, as the DefenseFinder wrapper already does. That gives ordered, per-locus calls that can be intersected with islands and plasmids.

## Fine as is
- IslandPath-DIMOB, geNomad (plasmid/virus/provirus), CheckV and IntegronFinder outputs are usable.
- They are used in `hgt/` and `islands/`.
