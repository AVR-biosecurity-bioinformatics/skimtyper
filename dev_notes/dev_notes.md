

## Create Qfly test datasets
This test data set uses a small segment of Qfly chromosome 1: CM028320.1:50000-99999

Only read pairs where at least 1 of the reads aligns to this region are included. 
```
module load SAMtools/1.21-GCC-13.3.0
module load BEDTools/2.31.1-GCC-13.3.0

# Sample 1 EM6.bam
samtools view -b --fetch-pairs /group/pathogens/IAWS/Projects/Tephritid/Skim/bams_recal_nind1maf1/bams/EM6.bam "CM028320.1:50000-99999" "CM028321.1:50000-59999" \
| samtools sort -n > subset.bam
bedtools bamtofastq -i subset.bam -fq test_data/qfly/EM6_subset_R1.fastq -fq2 test_data/qfly/EM6_subset_R2.fastq
gzip -f test_data/qfly/EM6_subset_R1.fastq test_data/qfly/EM6_subset_R2.fastq

# Sample 2 EM3.bam
samtools view -b --fetch-pairs /group/pathogens/IAWS/Projects/Tephritid/Skim/bams_recal_nind1maf1/bams/EM3.bam "CM028320.1:50000-99999" "CM028321.1:50000-59999" \
| samtools sort -n > subset.bam
bedtools bamtofastq -i subset.bam -fq test_data/qfly/EM3_subset_R1.fastq -fq2 test_data/qfly/EM3_subset_R2.fastq
gzip -f test_data/qfly/EM3_subset_R1.fastq test_data/qfly/EM3_subset_R2.fastq

# Sample 3 F3.bam
samtools view -b --fetch-pairs /group/pathogens/IAWS/Projects/Tephritid/Skim/bams_recal_nind1maf1/bams/F3.bam "CM028320.1:50000-99999" "CM028321.1:50000-59999" \
 | samtools sort -n > subset.bam
bedtools bamtofastq -i subset.bam -fq test_data/qfly/F3_subset_R1.fastq -fq2 test_data/qfly/F3_subset_R2.fastq
gzip -f test_data/qfly/F3_subset_R1.fastq test_data/qfly/F3_subset_R2.fastq

# Sample 4 F2xM12-F1.bam
samtools view -b --fetch-pairs /group/pathogens/IAWS/Projects/Tephritid/Skim/bams_recal_nind1maf1/bams/F2xM12-F1.bam "CM028320.1:50000-99999" "CM028321.1:50000-59999" \
 | samtools sort -n > subset.bam
bedtools bamtofastq -i subset.bam -fq test_data/qfly/F2xM12-F1_subset_R1.fastq -fq2 test_data/qfly/F2xM12-F1_subset_R2.fastq
gzip -f test_data/qfly/F2xM12-F1_subset_R1.fastq test_data/qfly/F2xM12-F1_subset_R2.fastq

# Subset reference genome to that portion - Fix header with sed to avoid error with gatk
samtools faidx /group/referencedata/mspd-db/genomes/insect/bactrocera_tryoni/GCA_016617805.2_CSIRO_BtryS06_freeze2_genomic.fna "CM028320.1:50000-99999" "CM028321.1:50000-59999" | sed 's/:.*$//g' > test_data/qfly/test_qfly_genome.fa

# add mitochondrial genome to reference genome
cat /group/referencedata/mspd-db/genomes/insect/bactrocera_tryoni/mitogenome/HQ130030.1_Bactrocera_tryoni_mitochondrion.fa >> test_data/qfly/test_qfly_genome.fa

# create sample data sheet
fwd=$( find test_data/qfly/ -maxdepth 1 -name '*.fastq.gz' -type f | grep '_R1' | sort | uniq )
rev=$(echo "$fwd" | sed 's/_R1/_R2/g' )
sample_id=$(echo "$fwd" | sed 's/_subset.*$//g' | sed 's/^.*\///g')

# Create fake population labels
pop=$(echo -e "Pop1\nPop1\nPop2\nPop3")

# format sample,fastq_1,fastq_2,
paste -d ',' <(echo "sample_id") <(echo "pop") <(echo "fwd") <(echo "rev") > test_data/qfly/test_samplesheet.csv
paste -d ',' <(echo "$sample_id") <(echo "$pop")  <(echo "$fwd") <(echo "$rev") >> test_data/qfly/test_samplesheet.csv
```

## Create reference panel

NOTE: Panel coordinates need to be lifted over to subset reference genome
```
INVCF="/group/pathogens/IAWS/Personal/Alexp/skimseq_qfly/output/results/vcf/filtered_snp/snp.vcf.gz"
REGION="CM028320.1:50000-99999"
CONTIG="CM028320.1"
START=50000
OUT_PREFIX="qfly_panel.CM028320.1_50000_99999.filtered"

SUBSET_REF="test_data/qfly/test_qfly_genome.fa"

# 1) Extract region and drop explicitly named samples
bcftools view \
  -r "${REGION}" \
  -s ^EM3,EM6,F3,F2xM12-F1 \
  -Oz \
  -o subset.tmp.vcf.gz \
  "${INVCF}"

bcftools index -t subset.tmp.vcf.gz

# 2) Compute per-sample missingness on this subset
bcftools query -f'[%SAMPLE\t%GT\n]' subset.tmp.vcf.gz \
| awk '
BEGIN { OFS="\t" }
{
    total[$1]++
    if ($2 ~ /\./) missing[$1]++
}
END {
    for (s in total) {
        miss = (s in missing ? missing[s] : 0)
        frac = miss / total[s]
        print s, frac, miss, total[s]
    }
}' \
| sort -k2,2nr > sample_missingness.tsv

# 3) Drop samples with >90% missing in this subset
awk '$2 > 0.9 {print $1}' sample_missingness.tsv > drop.missing90.samples.txt

# 4) Build minimal filtered panel in ORIGINAL genome coordinates first
bcftools view \
  -S ^drop.missing90.samples.txt \
  -Ou \
  subset.tmp.vcf.gz \
| bcftools annotate \
  -Ou \
  -x INFO,^FORMAT/GT \
| bcftools +fill-tags \
  -Oz \
  -o "${OUT_PREFIX}.origcoords.vcf.gz" \
  -- -t AC,AN,AF,NS

bcftools index -t "${OUT_PREFIX}.origcoords.vcf.gz"

# 5) Lift panel positions to SUBSET-reference coordinates
#    new_POS = old_POS - START + 1
bcftools view "${OUT_PREFIX}.origcoords.vcf.gz" \
| awk -v contig="${CONTIG}" -v start="${START}" 'BEGIN{OFS="\t"}
/^##/ { print; next }
/^#CHROM/ { print; next }
{
    if ($1 == contig) {
        $2 = $2 - start + 1
    }
    print
}' \
| bgzip -c > "${OUT_PREFIX}.subsetcoords.vcf.gz"

bcftools index -t "${OUT_PREFIX}.subsetcoords.vcf.gz"

# 6) Update contig lengths/header to match the subset FASTA
#    (only if you will call/merge against the subset reference)
samtools faidx "${SUBSET_REF}"

bcftools reheader \
  -f "${SUBSET_REF}.fai" \
  -o panel.vcf.gz \
  "${OUT_PREFIX}.subsetcoords.vcf.gz"

bcftools index -t panel.vcf.gz


```

### Run test datasets


Run the Qfly test dataset using the test profile
```
module purge
module load Java/17

export NXF_VER=23.05.0-edge

# Run tests on local node 
nextflow run . -profile debug,test -resume

```


# Test aquilonis

```
cd /group/pathogens/IAWS/Personal/Alexp/skimtyper

ml BCFtools/1.23.1-GCC-13.3.0

awk 'NR>1 {print $1}' popmap_ref.tsv > samples_to_keep.txt


bcftools view \
    --threads 8 \
    -S samples_to_keep.txt \
    -Ou \
  /group/pathogens/IAWS/Personal/Alexp/skimseq_qfly/vcf_in_sitelist.vcf.gz \
  | bcftools annotate \
  -Ou \
  -x INFO,^FORMAT/GT \
| bcftools +fill-tags \
  -Oz \
  -o qfly_ref_panel.vcf.gz \
  -- -t AC,AN,AF,NS

bcftools index -t qfly_ref_panel.vcf.gz



# One or more FASTQ roots:
fqroots=(
  /group/pathogens/IAWS/Projects/Tephritid/Skim/fastq
  /group/sequencing/MGI/250723_V350314028/Project_PATHOGENS
  /group/sequencing/MGI/V350246164/Project_PATHOGENS
  /group/sequencing/MGI/V350246099/L04
  /group/sequencing/MGI/251127_R2130410240006_V350358133/Project_PATHOGENS
  /group/sequencing/MGI/250918_R2130410240006_V350353273/Project_PATHOGENS
  /group/sequencing/MGI/260129_R2130410240006_V350384689/Project_PATHOGENS
  )
  
out=sample_sheet.csv
printf "sample,pop,fwd,rev\n" > sample_sheet.csv

# read popmap (skip header)
tail -n +2 popmap_perth.tsv | while IFS=$'\t' read -r sample pop; do
  # Find all matching R1 files across ALL roots (NUL-delimited for safety)
  while IFS= read -r -d '' fwd; do
    rev="${fwd/_R1_/_R2_}"
    if [[ ! -f "$rev" ]]; then
      printf 'WARN: missing R2 for %s (%s)\n' "$sample" "$fwd" >&2
      continue
    fi
    printf "%s,%s,%s,%s\n" "$sample" "$pop" "$fwd" "$rev" >> "$out"
  done < <(
    # Adjust the name patterns to your naming scheme
    find "${fqroots[@]}" -maxdepth 3 -type f -name '*.fastq.gz' \
         -name '*_R1_*' -name "*${sample}_*" -print0 \
    | sort -z -u
  )
done
```