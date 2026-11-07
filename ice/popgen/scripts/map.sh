#!/bin/bash
# map.sh SAMPLE R1 R2 NPAIRS : first NPAIRS read pairs -> bwa mem vs all 911 MAGs -> MAPQ>=20, primary, >=95% identity, >=80% of read aligned -> sorted BAM
set -euo pipefail
S=$1; R1=$2; R2=$3; N=${4:-2000000}
D=/data/dev/rnaseq-exploration/ice/popgen
OUT=$D/bam/$S.bam
[ -s $OUT ] && exit 0
bwa mem -t 14 -v 1 $D/ref/allbins.fa <(zcat $R1 | head -n $((4*N))) <(zcat $R2 | head -n $((4*N))) 2> $D/bam/$S.bwa.log \
 | samtools view -h -q 20 -F 0x904 - \
 | awk 'BEGIN{OFS="\t"} /^@/{print;next} {nm=-1; for(i=12;i<=NF;i++) if($i~/^NM:i:/){nm=substr($i,6)+0}
     c=$6; al=0; ql=0; while(match(c,/[0-9]+[MIDSH=X]/)){n=substr(c,RSTART,RLENGTH-1)+0; o=substr(c,RSTART+RLENGTH-1,1); if(o~/[M=XD]/)al+=n; if(o~/[M=XI]/)ql+=n; if(o~/[MIS=X]/)rl+=0; c=substr(c,RSTART+RLENGTH)}
     if(nm<0||al==0) next; if((al-nm)/al<0.95) next; if(ql<0.8*length($10)) next; print}' \
 | samtools sort -@4 -m 1G -o $OUT.tmp - && mv $OUT.tmp $OUT && samtools index $OUT
samtools idxstats $OUT | awk -v s=$S '{m+=$3} END{print s"\t"m}' >> $D/bam/mapped_counts.tsv
