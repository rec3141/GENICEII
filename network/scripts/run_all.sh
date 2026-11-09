#!/bin/bash
# Full network rebuild, in order. Rerun after bins, EE MAGs (data/hull_mags.tsv), organelle groups or vOTUs change.
# Julia/FlashWeave: conda env /data/scratch/grid-envs/conda/julia, depot /data/scratch/grid-envs/julia-depot.
# R: SpiecEasi (system library); validate_spring.py also needs latentcor in /data/scratch/grid-envs/Rlib and a SPRING checkout.
set -euo pipefail
cd "$(dirname "$0")/.."
PY=/data/miniforge3/bin/python3
$PY scripts/nodes.py                      # nodes, contig assignment, DNA/RNA tables (cache/)
$PY scripts/export.py                     # cache/{dna,rna}_X.tsv, *_meta.tsv
$PY scripts/run_spearman.py               # design-residual Spearman
$PY -u scripts/run_spring.py              # SPRING (Python port)
Rscript scripts/run_spieceasi.R .         # SPIEC-EASI mb
JULIA_DEPOT_PATH=/data/scratch/grid-envs/julia-depot /data/scratch/grid-envs/conda/julia/bin/julia -p 6 scripts/run_flashweave.jl .
$PY scripts/consensus.py                  # consensus, overlaps, known pairs
$PY scripts/layout.py                     # t-SNE positions
$PY ../scripts/lineage.py                 # site/lineage.json (adds the vOTU nodes)
$PY scripts/network_site.py               # site/network.json + network_{dna,rna}.json (full edge sets, permissivity steps)
# validation of the SPRING port (simulation + R reference): $PY scripts/validate_spring.py
