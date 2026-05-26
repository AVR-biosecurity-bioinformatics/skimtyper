#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = ref_genome
# $4 = sample
# $5 = interval hash
# $6 = targets_file

## Parse positional input args, the rest are xported
CPUS="${1}"
MEM_GB="${2}"
REF="${3}"
SAMPLE="${5}"
IHASH="${5}"
PANEL_VCF="${6}"

# set up bcftools filter flags
if [[ "${RMDUP}" == "false" ]]; then
  FILTER_FLAGS="--ns DUP -G UNMAP,SECONDARY,QCFAIL"
else
  FILTER_FLAGS="--ns UNMAP,SECONDARY,QCFAIL,DUP"
fi

# Call variants on target sites only with mpleup
bcftools mpileup \
    --threads ${CPUS} \
    --bam-list cram.list \
    --max-depth ${MAXDEPTH} \
    --fasta-ref ${REF} \
    --min-BQ ${MINBQ} \
    --min-MQ ${MINMQ} \
    --regions-file ${PANEL_VCF} \
    ${FILTER_FLAGS} \
    --annotate FORMAT/DP,FORMAT/AD,INFO/AD \
    --indels-cns \
    --indel-size 110 \
    | bcftools call \
    -Ou \
    -a FORMAT/GP,FORMAT/GQ \
    --ploidy ${PLOIDY} \
    --constrain alleles \
    --regions-file ${PANEL_VCF} \
    --insert-missed \
    --multiallelic-caller \
    --prior ${MUTATION_RATE} \
  | bcftools +setGT \
    -Ou -- \
    -t q -n . -i 'FMT/DP=0' \
  | bcftools annotate \
    --threads "${CPUS}" \
    --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' \
    -Oz9 -o "${SAMPLE}.${IHASH}.vcf.gz"

# index output
bcftools index -t ${SAMPLE}.${IHASH}.vcf.gz

# Extra annotatiomns that can be added
# mpileup: FORMAT/SP
# call: FORMAT/GP,INFO/PV4
