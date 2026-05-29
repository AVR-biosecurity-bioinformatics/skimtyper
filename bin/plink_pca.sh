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
    --mem ${2} \
    --bfile ${outname} \
    --snps-only just-acgt \
    --max-alleles 2 \
    --mac ${MIN_MAC} \
    --maf ${MIN_MAF} \
    --geno ${SITE_MAX_MISSING} \
    --mind ${SAMPLE_MAX_MISSING} \
    --allow-extra-chr \
    --make-bed \
    --out ${outname}.filtered

# Then LD prune
plink2 \
    --threads ${1} \
    --mem ${2} \
    --bfile ${outname}.filtered \
    --allow-extra-chr \
    --indep-pairwise ${LD_WINSIZE} ${LD_STEPSIZE} ${LD_THRESHOLD} \
    --out ${outname}.pruned

# Then run PCA on filtered and pruned snps
plink2 \
    --threads ${1} \
    --mem ${2} \
    --bfile ${outname}.filtered \
    --extract ${outname}.filtered.prune.in \
    --allow-extra-chr \
    --pca 10 \
    --out ${outname}.pca
