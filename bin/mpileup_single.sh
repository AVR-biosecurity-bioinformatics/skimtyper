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
# $7 = cram

## Parse positional input args, the rest are xported
CPUS="${1}"
MEM_GB="${2}"
REF="${3}"
SAMPLE="${4}"
IHASH="${5}"
PANEL_VCF="${6}"
CRAM="${7}"

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

## Buid annotation table of AN and AC for ref panel priors
bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t%AN\t%AC\t%ID\n' "${PANEL_VCF}" \
  | bgzip -c > panel.acan.tsv.gz
tabix -f -s1 -b2 -e2 panel.acan.tsv.gz

cat > panel.acan.hdr <<'EOF'
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes from reference panel">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Alternate allele count from reference panel">
##INFO=<ID=PANEL_ID,Number=1,Type=String,Description="Reference panel SNP ID">
EOF

# flag to use panel priors or not
if [[ "${REF_PANEL_PRIORS}" == "false" ]]; then
  PRIOR_FLAGS="--prior ${MUTATION_RATE}"
else
  PRIOR_FLAGS="--prior-freqs AN,AC"
fi

# Call variants on target sites only with mpleup
bcftools mpileup \
    -Ou \
    --threads "${CPUS}" \
    --max-depth "${MAXDEPTH}" \
    --fasta-ref "${REF}" \
    --min-BQ "${MINBQ}" \
    --min-MQ "${MINMQ}" \
    --regions-file "${PANEL_VCF}" \
    ${FILTER_FLAGS} \
    --annotate FORMAT/DP,FORMAT/AD,INFO/AD \
    --indels-cns \
    --indel-size 110 \
    "${CRAM}" \
| bcftools call \
    -Ou \
    -a FORMAT/GP,FORMAT/GQ \
    -A \
    --ploidy "${PLOIDY}" \
    --constrain alleles \
    -T panel.alleles.tsv.gz \
    --insert-missed \
    --multiallelic-caller \
    ${PRIOR_FLAGS} \
| bcftools +setGT \
    -Ou -- \
    -t q -n . -i 'FMT/DP=0' \
| bcftools annotate \
    -a panel.acan.tsv.gz \
    -h panel.acan.hdr \
    -c CHROM,POS,REF,ALT,INFO/AN,INFO/AC,ID \
    -Oz9 -o "${SAMPLE}.${IHASH}.vcf.gz"

# index output
bcftools index -t ${SAMPLE}.${IHASH}.vcf.gz

# Extra annotatiomns that can be added
# mpileup: FORMAT/SP
# call: FORMAT/GP,INFO/PV4
