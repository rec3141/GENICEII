#!/bin/bash
#SBATCH --account=def-rec3141 --time=1:00:00 --mem=8G --cpus-per-task=1 --job-name=islands --output=islands-%j.log
cd ~/scratch/project_rnaseq/explore/islands
B=~/scratch/project_nanopore/mag_marine/annotation/bakta/extra/annotation.tsv
cut -f2 genomic_islands.tsv | tail -n +2 | sort -u > island_contigs.txt
awk -F"\t" "NR==FNR{k[\$1]=1;next} (\$1 in k)" island_contigs.txt $B | gzip > island_contig_features.tsv.gz
SIF=/home/rec3141/GENICE/danaSeq/mag_analysis/.danaseq-mag-analysis.sif
/opt/singularity-ce-4.4.1/bin/apptainer exec -B $(readlink -f ~/scratch) $SIF bash -c "source /opt/conda/etc/profile.d/conda.sh && conda activate /opt/conda/envs/dana-mag-assembly && python3 gc.py $HOME/scratch/project_nanopore/assembly_marine/assembly/assembly.fasta"
echo done
