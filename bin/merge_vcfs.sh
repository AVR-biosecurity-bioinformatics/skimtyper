#!/bin/bash
set -euo pipefail
## args:
# $1 = cpus
# $2 = mem (GB)
# $3 = ref_genome
# $4 = outname

## Parse positional input args
CPUS="${1}"
MEM_GB="${2}"
REF="${3}"
OUTNAME="${4}"

# Detect vcf type (.g.vcf.gz or .vcf.gz) from the first file
first=$(head -n1 vcf.list || true)
if [[ -z "$first" ]]; then
    echo "No VCFs found (expected *.vcf.gz)" >&2
    exit 1
fi

# Exit if extension not recognised
if [[ "$first" == *.g.vcf.gz ]]; then
    ext=".g.vcf.gz"
    GVCF_FLAG="--gvcf ${REF}"
elif [[ "$first" == *.vcf.gz ]]; then
    ext=".vcf.gz"
    GVCF_FLAG=""
else
    echo "File extension not recognised: $first" >&2
    exit 1
fi

bcftools merge \
  -Ou \
  --threads ${CPUS} \
  ${GVCF_FLAG} \
  --file-list vcf.list \
  --force-samples \
  | bcftools +fill-tags \
  -Oz \
  -o ${OUTNAME}.merged.${ext} \
  -- -t AC,AN,AF,NS

# Index output
bcftools index -t --threads ${CPUS} ${OUTNAME}.merged.${ext}
