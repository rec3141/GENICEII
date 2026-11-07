#!/bin/bash
# Map every sample in samples_map.tsv (first 2M pairs) against the competitive all-MAG index; per-sample wall time -> bam/times.log
D=/data/dev/rnaseq-exploration/ice/popgen
[ -s $D/ref/allbins.fa.sa ] || { echo "index missing"; exit 1; }
tail -n +2 $D/samples_map.tsv | awk -F'\t' '{print $1"\t"$7"\t"$9}' | while IFS=$'\t' read s r1 r2; do
  t0=$SECONDS
  if $D/scripts/map.sh "$s" "$r1" "$r2" 2000000; then echo "$s $((SECONDS - t0)) s" >> $D/bam/times.log; else echo "FAIL $s" >> $D/bam/times.log; fi
done
echo done
