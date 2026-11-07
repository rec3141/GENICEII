#!/bin/bash
#SBATCH --account=def-rec3141 --time=3:00:00 --mem=32G --cpus-per-task=16 --job-name=viruses --output=viruses-%j.log
cd ~/scratch/project_rnaseq/explore/viruses
SIF=/home/rec3141/GENICE/danaSeq/mag_analysis/.danaseq-mag-analysis.sif
M=~/scratch/project_nanopore
AP="/opt/singularity-ce-4.4.1/bin/apptainer exec -B $(readlink -f ~/scratch) $SIF bash -c"
$AP "source /opt/conda/etc/profile.d/conda.sh && conda activate /opt/conda/envs/dana-mag-assembly && samtools faidx $M/assembly_marine/assembly/assembly.fasta -r viral_ids.txt > viral.fna && skani triangle --sparse -t 16 -m 200 viral.fna > skani_sparse.tsv 2> skani.log"
echo skani done
$AP "source /opt/conda/etc/profile.d/conda.sh && conda activate /opt/conda/envs/dana-mag-assembly && python3 vdepth.py $M/mag_marine/mapping/depths.txt"
echo done
