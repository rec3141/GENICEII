#!/bin/bash
#SBATCH --account=def-rec3141 --time=2:00:00 --mem=32G --cpus-per-task=16 --job-name=vseq --output=vseq-%j.log
cd ~/scratch/project_rnaseq/explore/viruses
SIF=/home/rec3141/GENICE/danaSeq/mag_analysis/.danaseq-mag-analysis.sif
M=~/scratch/project_nanopore
/opt/singularity-ce-4.4.1/bin/apptainer exec -B $(readlink -f ~/scratch) $SIF bash -c "source /opt/conda/etc/profile.d/conda.sh && conda activate /opt/conda/envs/dana-mag-assembly && skani triangle -i --sparse -t 16 -m 200 viral_all.fna > skani_sparse.tsv 2> skani.log && wc -l skani_sparse.tsv"
echo done
