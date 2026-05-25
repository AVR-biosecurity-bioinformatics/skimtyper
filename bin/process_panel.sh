#!/bin/bash
set -euo pipefail
## args:
# $1 = cpus
# $2 = mem (GB)
# $3 = vcf
# $4 = outname

bcftools view --threads ${1} -G -O z9 -o panel.sites.vcf.gz ${3}
bcftools index -t --threads ${1} "panel.sites.vcf.gz"