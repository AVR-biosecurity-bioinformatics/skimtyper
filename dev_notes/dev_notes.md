## Create reference panel
```
bcftools view \
  -r CM028320.1:50000-99999 \
  -s ^EM3,EM6,F3,F2xM12-F1 \
  -Oz \
  -o subset.tmp.vcf.gz \
  /group/pathogens/IAWS/Personal/Alexp/skimseq_qfly/output/results/vcf/filtered_snp/snp.vcf.gz

bcftools index -t subset.tmp.vcf.gz


bcftools query -f'[%SAMPLE\t%GT\n]' subset.tmp.vcf.gz \
| awk '
BEGIN { OFS="\t" }
{
    total[$1]++
    if ($2=="./." || $2==".|.") missing[$1]++
}
END {
    for (s in total) {
        miss = (s in missing ? missing[s] : 0)
        frac = miss / total[s]
        print s, frac, miss, total[s]
    }
}' \
| sort -k2,2nr > sample_missingness.tsv

awk '$2 > 0.9 {print $1}' sample_missingness.tsv > drop.missing90.samples.txt


bcftools view \
  -S ^drop.missing90.samples.txt \
  -Ou \
  subset.tmp.vcf.gz \
| bcftools annotate \
  -Ou \
  -x INFO,^FORMAT/GT \
| bcftools +fill-tags \
  -Oz \
  -o qfly_panel.CM028320.1_50000_99999.filtered.vcf.gz \
  -- -t AC,AN,AF,NS

bcftools index -t qfly_panel.CM028320.1_50000_99999.filtered.vcf.gz

```

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

### Run test datasets


Run the Qfly test dataset using the test profile
```
module purge
export NXF_VER=23.05.0-edge
module load Java/17

# Run tests on local node 
nextflow run . -profile debug,test -resume

```