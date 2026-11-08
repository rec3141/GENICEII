#!/bin/bash
# One-shot fetch of existing danaSeq mag_marine evidence for identify.py (2026-11-08). ONE ssh call, no loop, no retry.
# Copies (read-only source): eukaryotic/marferret/*, eukaryotic/metaeuk/metaeuk_headers.tsv, taxonomy/rrna/rrna_{genes,contigs}.tsv,
# rrna_sequences.fasta; and the sequences of the contigs listed in $1 (samtools faidx in the danaSeq container) as contigs.fa.
# Usage: fetch_grex.sh contig_list.txt  -> out/identify/grex/
set -euo pipefail
OUT=/data/dev/rnaseq-exploration/euks/mags/out/identify/grex; mkdir -p "$OUT"
REMOTE='set -e; export PATH=/opt/singularity-ce-4.4.1/bin:$PATH; M=~/scratch/project_nanopore/mag_marine
R=~/scratch/project_rnaseq/marine/references/marine_coassembly; SIF=~/GENICE/danaSeq/mag_analysis/.danaseq-mag-analysis.sif
T=$(mktemp -d ~/scratch/identify_fetch.XXXX); B=$(readlink -f ~/scratch); cat > $T/list.txt
apptainer exec -B $B $SIF bash -c "source /opt/conda/etc/profile.d/conda.sh && conda activate /opt/conda/envs/dana-mag-assembly && samtools faidx $R/marine_coassembly.fasta -r $T/list.txt > $T/contigs.fa" >&2
cp $M/eukaryotic/marferret/marferret_proteins.tsv $M/eukaryotic/marferret/marferret_contigs.tsv $M/eukaryotic/metaeuk/metaeuk_headers.tsv \
   $M/taxonomy/rrna/rrna_genes.tsv $M/taxonomy/rrna/rrna_contigs.tsv $M/taxonomy/rrna/rrna_sequences.fasta $T/
tar -C $T -czf - . ; rm -rf $T'
ssh -o BatchMode=yes concentration "ssh -o BatchMode=yes grex $(printf '%q' "bash -c $(printf '%q' "$REMOTE")")" < "$1" | tar -C "$OUT" -xzf -
ls -la "$OUT"
