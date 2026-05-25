#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = mem
# $3 = vcf
# $4 = counts_per_chunk

TARGET_COUNTS_PER_CHUNK=$(awk -v x="${4}" 'BEGIN {printf("%d\n",x)}')
bcftools +scatter -n $TARGET_COUNTS_PER_CHUNK ${3} -o . -p chunk_

# Rename each output file
for i in *chunk_*.vcf;do
  # Pad output chunk names
  n=$(basename "$i" | sed -E 's/^chunk_([0-9]+)\.vcf/\1/')
  pad=$(printf "%05d" "$n")

  # compute hash of this chunk’s contents
  hash=$(md5sum "$i" | awk '{print $1}')

  # Output just column 1:4 as a gzipped bed
  out="${pad}${hash}.vcf.gz"
  cat "$i" | bgzip -c > "$out"
  bcftools index -t "$out"

  rm $i
done
