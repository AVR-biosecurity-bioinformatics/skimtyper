/*
    Genotype samples using GATK
*/

//// import modules
include { MERGE_VCFS                                             } from '../modules/merge_vcfs' 
include { CREATE_INTERVAL_CHUNKS as CREATE_INTERVAL_CHUNKS_MP    } from '../modules/create_interval_chunks'
include { MPILEUP                                                } from '../modules/mpileup'

workflow MPILEUP_CALLING {

    take:
    ch_sample_names
    ch_sample_cram
    ch_genome_indexed
    ch_sites_to_genotype
    ch_read_counts

    main: 

    // combine sample-level cram with each interval_bed file and interval chunk
    ch_sample_cram 
        .combine ( ch_sites_to_genotype )
        .map { sample, cram, crai, interval_hash, interval_bed, bed_tbi, sites_vcf, sites_tbi -> [ interval_hash, cram, crai ] }
        .groupTuple ( by: [0,1] )
        // join to get back interval_file
        .join ( ch_sites_to_genotype, by: [0,1] )
        // variant type and interval hash columns are combined into a single string for compatibility with mpileup
        .map { interval_hash, cram, crai, interval_bed, bed_tbi, sites_vcf, sites_tbi -> tuple(interval_hash, sites_vcf, sites_tbi, cram, crai) }
	    .set { ch_cram_to_genotype }

    // Calculate cohort size for memory scaling
    ch_cohort_size = ch_sample_names.unique().count()

    // Call just target sites using mpileup
    MPILEUP (
        ch_cram_to_genotype,
        ch_genome_indexed,
        ch_cohort_size
    )
     
    // Merge seperate VCFs
    MPILEUP.out.vcf
        .map { interval_chunk, interval_bed, bed_tbi, vcf, tbi -> tuple('unfiltered', vcf, tbi) }
        .groupTuple(by: 0)
        .set { ch_vcf_to_merge }

    MERGE_VCFS (
        ch_vcf_to_merge
    )

    emit: 
    vcf = MPILEUP.out.vcf

}