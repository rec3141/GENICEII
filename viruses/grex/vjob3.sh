#!/bin/bash
#SBATCH --account=def-rec3141 --time=1:00:00 --mem=8G --cpus-per-task=1 --job-name=vfeat --output=vfeat-%j.log
cd ~/scratch/project_rnaseq/explore/viruses
E=~/scratch/project_nanopore/mag_marine/annotation/bakta/extra
awk -F'\t' 'NR==FNR{k[$1]=1;next} ($1 in k)' provirus_contigs.txt $E/annotation.tsv > provirus_features.tsv
cut -f6 provirus_features.tsv | grep . | sort -u > tags.txt
awk -F'\t' 'NR==FNR{k[$1]=1;next} ($2 in k)' tags.txt $E/locus_tag_map.tsv > provirus_tagmap.tsv
grep -h "^contig\|^#" $E/../../../mge/checkv/quality_summary.tsv | head -0
wc -l provirus_features.tsv provirus_tagmap.tsv
echo done
