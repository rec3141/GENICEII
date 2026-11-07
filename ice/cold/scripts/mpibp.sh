#!/bin/bash
#SBATCH --account=def-rec3141 --time=1:00:00 --mem=8G --cpus-per-task=4 --job-name=mpibp --output=%x-%j.log
set -euo pipefail
export PATH=/opt/singularity-ce-4.4.1/bin:$PATH
SIF=/home/rec3141/GENICE/danaSeq/mag_analysis/.danaseq-mag-analysis.sif
SC=$(readlink -f ~/scratch); O=$SC/project_rnaseq/explore/ice/cold; M=$SC/project_nanopore/mag_marine
cd $O
awk 'NR==FNR{k[">"$1]=1;next} /^>/{p=($1 in k)} p' marinomonas_ids.txt $M/annotation/bakta/extra/annotation.faa > marinomonas.faa
grep -c ">" marinomonas.faa
apptainer exec --bind $SC $SIF bash -c "source /opt/conda/etc/profile.d/conda.sh && conda activate /opt/conda/envs/dana-mag-annotate && diamond makedb --in marinomonas.faa -d mar -p 4 --quiet && diamond blastp -q mpibp.faa -d mar -o mpibp.m8 -p 4 --more-sensitive -e 1e-5 --max-target-seqs 50 --outfmt 6 qseqid sseqid pident length qstart qend sstart send evalue bitscore slen --quiet"
cat mpibp.m8
