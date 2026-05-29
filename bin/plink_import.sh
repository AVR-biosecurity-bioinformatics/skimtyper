#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = mem
# $3 = outname
# $4 = vcf

# Create PLINK bed file
plink2 \
    --threads ${1} \
    --memory ${2} \
    --vcf ${4} \
    --allow-extra-chr \
    --double-id \
    --make-bed \
    --out ${3}