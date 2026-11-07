#!/bin/bash
#SBATCH --account=def-rec3141 --time=6:00:00 --mem=48G --cpus-per-task=32 --job-name=icecold --output=%x-%j.log
set -euo pipefail
export PATH=/opt/singularity-ce-4.4.1/bin:$PATH
SIF=/home/rec3141/GENICE/danaSeq/mag_analysis/.danaseq-mag-analysis.sif
SC=$(readlink -f ~/scratch); O=$SC/project_rnaseq/explore/ice/cold; M=$SC/project_nanopore/mag_marine
P=$SC/databases/eggnog_db/pfam/Pfam-A.hmm
A="source /opt/conda/etc/profile.d/conda.sh && conda activate /opt/conda/envs/dana-mag-annotate"
cd $O
apptainer exec --bind $SC $SIF bash -c "$A && hmmfetch -f $P pfam_list.txt > cold.hmm && hmmsearch --cpu 32 --cut_ga --noali --tblout pf.tbl --domtblout pf.domtbl -o /dev/null cold.hmm $M/annotation/bakta/extra/annotation.faa"
apptainer exec --bind $SC $SIF bash -c "$A && python3 aacount.py $M/annotation/bakta/extra/annotation.faa $O && python3 bingc.py $M/binning/dastool/bins > bin_gc.tsv"
ls -la $O
