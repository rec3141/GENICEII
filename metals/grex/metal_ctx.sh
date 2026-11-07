#!/bin/bash
#SBATCH --account=def-rec3141 --time=1:00:00 --mem=16G --cpus-per-task=1 --job-name=metalctx --output=metalctx-%j.log
# Features on contigs that carry metal/toxicity genes (extra Bakta run, reconciled tags) + old->new tag map for those genes.
# Inputs (uploaded): contigs.txt, old_tags.txt (metal genes + confident HC genes, pre-Nov-2 extra tags).
cd ~/scratch/project_rnaseq/explore/metals
B=~/scratch/project_nanopore/mag_marine/annotation/bakta/extra
awk -F'\t' 'NR==FNR{k[$1]=1;next} ($1 in k)' old_tags.txt $B/locus_tag_map.tsv > tagmap_subset.tsv
awk -F'\t' 'NR==FNR{k[$1]=1;next} ($1 in k)' contigs.txt $B/annotation.tsv | cut -f1-8 | gzip > contig_features.tsv.gz
wc -l tagmap_subset.tsv old_tags.txt; zcat contig_features.tsv.gz | wc -l
echo done
