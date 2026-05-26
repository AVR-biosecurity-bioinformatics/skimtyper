#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = ref_genome
# $4 = sample
# $5 = interval hash
# $6 = panel_vcf

## Parse positional input args, the rest are xported
CPUS="${1}"
MEM_GB="${2}"
REF="${3}"
SAMPLE="${4}"
IHASH="${5}"
PANEL_VCF="${6}"

# set up bcftools filter flags
if [[ "${RMDUP}" == "false" ]]; then
  FILTER_FLAGS="--ns DUP -G UNMAP,SECONDARY,QCFAIL"
else
  FILTER_FLAGS="--ns UNMAP,SECONDARY,QCFAIL,DUP"
fi

 ## Build allele targets: CHROM POS REF,ALT  (tabix indexed)
bcftools view -m2 -M2 "${PANEL_VCF}" \
    | bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' \
    | bgzip -c > panel.alleles.tsv.gz
tabix -s1 -b2 -e2 panel.alleles.tsv.gz

## Buid ref panel priors if set
if [[ "${REF_PANEL_PRIORS:-false}" == "true" ]]; then
  bcftools view -m2 -M2 "${PANEL_VCF}" \
    | bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t%AN\t%AC\n' \
    | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,($5-$6) "," $6}' \
    | bgzip -c > panel.prior_freqs.tsv.gz

  tabix -f -s1 -b2 -e2 panel.prior_freqs.tsv.gz
  PRIOR_ARG=(--prior-freqs panel.prior_freqs.tsv.gz)
else
  PRIOR_ARG=(--prior "${MUTATION_RATE}")
fi

# Call variants on target sites only with mpleup
bcftools mpileup \
    -Ou \
    --threads ${CPUS} \
    --max-depth ${MAXDEPTH} \
    --fasta-ref ${REF} \
    --min-BQ ${MINBQ} \
    --min-MQ ${MINMQ} \
    --regions-file ${PANEL_VCF} \
    ${FILTER_FLAGS} \
    --annotate FORMAT/DP,FORMAT/AD,INFO/AD \
    --indels-cns \
    --indel-size 110 \
    ${SAMPLE}.cram \
    | bcftools call \
    -Ou \
    -a FORMAT/GP,FORMAT/GQ \
    --ploidy ${PLOIDY} \
    --constrain alleles \
    -T panel.alleles.tsv.gz \
    --insert-missed \
    --multiallelic-caller \
    ${PRIOR_ARG} \
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
