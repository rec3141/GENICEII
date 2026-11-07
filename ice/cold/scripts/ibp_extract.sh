#!/bin/bash
#SBATCH --account=def-rec3141 --time=1:00:00 --mem=16G --cpus-per-task=4 --job-name=ibp-extract --output=%x-%j.log
set -euo pipefail
export PATH=/opt/singularity-ce-4.4.1/bin:$PATH
M=~/scratch/project_nanopore/mag_marine; R=~/scratch/project_rnaseq/marine/references/marine_coassembly
SIF=~/GENICE/danaSeq/mag_analysis/.danaseq-mag-analysis.sif; B=$(readlink -f ~/scratch)
cd ~/scratch/project_rnaseq/explore/ice/cold
awk "NR==FNR{k[\$1];next} /^>/{split(substr(\$1,2),a,\" \"); p=(a[1] in k)} p" ibp_gene_ids.txt $M/annotation/bakta/extra/annotation.faa > ibp_prot.faa
grep -c ">" ibp_prot.faa
zcat $R/marine_coassembly.genes.tsv.gz | awk -F"\t" "NR==FNR{c[\$1];next} FNR==1||(\$2 in c)" ibp_contigs.txt - | gzip > ibp_contig_genes.tsv.gz
awk -F"\t" "NR==FNR{c[\$1];next} FNR==1||(\$1 in c)" ibp_contigs.txt $M/assembly/gc.tsv > ibp_contig_gc.tsv
# +-20 kb windows around each IBP gene (clipped at contig ends by samtools)
zcat ibp_contig_genes.tsv.gz | awk -F"\t" "NR==FNR{k[\$1];next} (\$1 in k){s=\$3-20000; if(s<1)s=1; print \$2\":\"s\"-\"\$4+20000\"\t\"\$1}" ibp_gene_ids.txt - > ibp_regions.tsv
cut -f1 ibp_regions.tsv > ibp_regions.txt
apptainer exec -B $B $SIF bash -c "source /opt/conda/etc/profile.d/conda.sh && conda activate /opt/conda/envs/dana-mag-annotate && \
  mmseqs easy-search ibp_prot.faa ibp_prot.faa ibp_ava.m8 tmp_s --min-seq-id 0.3 -c 0.5 --cov-mode 0 --threads 4 --format-output query,target,fident,alnlen,qcov,tcov,evalue,bits >/dev/null"
apptainer exec -B $B $SIF bash -c "source /opt/conda/etc/profile.d/conda.sh && conda activate /opt/conda/envs/dana-mag-assembly && samtools faidx $R/marine_coassembly.fasta -r ibp_regions.txt > ibp_regions.fna"
ls $M/eukaryotic/marferret $M/eukaryotic/metaeuk 2>/dev/null | head
echo all-done
