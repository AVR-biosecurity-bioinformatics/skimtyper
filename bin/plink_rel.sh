#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = mem
# $3 = outname

# First apply filters
plink2 \
    --threads ${1} \
    --memory ${2} \
    --bfile ${3} \
    --snps-only just-acgt \
    --max-alleles 2 \
    --mac ${MIN_MAC} \
    --maf ${MIN_MAF} \
    --geno ${SITE_MAX_MISSING} \
    --mind ${SAMPLE_MAX_MISSING} \
    --allow-extra-chr \
    --make-bed \
    --out ${3}.filtered

# Then LD prune
plink2 \
    --threads ${1} \
    --memory ${2} \
    --bfile ${3}.filtered \
    --allow-extra-chr \
    --indep-pairwise ${LD_WINSIZE} ${LD_STEPSIZE} ${LD_THRESHOLD} \
    --out ${3}.filtered

# Then run PCA on filtered and pruned snps
plink2 \
    --threads ${1} \
    --memory ${2} \
    --bfile ${3}.filtered \
    --extract ${3}.filtered.prune.in \
    --allow-extra-chr \
    --make-rel square gz \
    --out ${3}
