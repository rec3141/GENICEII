#!/bin/bash
# Map the first 2M read pairs of each CMO pool metagenome (map/samples.tsv) to the vOTU representatives; keep primary alignments with
# MAPQ >= 1, >= 95% identity and >= 80% of the read aligned (as ice/popgen/scripts/map.sh); write reads per vOTU (no BAM kept).
D=/data/dev/rnaseq-exploration/viruses/map
tail -n +2 $D/samples.tsv | awk -F'\t' '{print $1"\t"$3"\t"$NF}' | while IFS=$'\t' read s r1 r2; do
  [ -s $D/counts/$s.tsv ] && continue
  bwa mem -t 14 -v 1 $D/votu_reps.fna <(zcat $r1 | head -n 8000000) <(zcat $r2 | head -n 8000000) 2> /dev/null \
   | samtools view -q 1 -F 0x904 - \
   | awk 'BEGIN{OFS="\t"} {nm=-1; for(i=12;i<=NF;i++) if($i~/^NM:i:/){nm=substr($i,6)+0}
       c=$6; al=0; ql=0; while(match(c,/[0-9]+[MIDSH=X]/)){n=substr(c,RSTART,RLENGTH-1)+0; o=substr(c,RSTART+RLENGTH-1,1); if(o~/[M=XD]/)al+=n; if(o~/[M=XI]/)ql+=n; c=substr(c,RSTART+RLENGTH)}
       if(nm<0||al==0) next; if((al-nm)/al<0.95) next; if(ql<0.8*length($10)) next; print $3}' \
   | sort | uniq -c | awk -v OFS='\t' '{print $2,$1}' > $D/counts/$s.tsv.tmp && mv $D/counts/$s.tsv.tmp $D/counts/$s.tsv
  echo "$s $(awk '{s+=$2} END{print s}' $D/counts/$s.tsv)" >> $D/map.log
done
echo ALLDONE >> $D/map.log
