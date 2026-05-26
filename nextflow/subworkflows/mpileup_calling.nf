/*
    Genotype samples using GATK
*/

//// import modules
include { MERGE_VCFS                                             } from '../modules/merge_vcfs' 
include { MPILEUP_SINGLE                                         } from '../modules/mpileup_single'
include { SCATTER_VCF                                            } from '../modules/scatter_vcf'

workflow MPILEUP_CALLING {

    take:
    ch_sample_names
    ch_sample_cram
    ch_genome_indexed
    ch_sites_to_genotype
    ch_read_counts

    main: 

    // TODO: Should i just do per-sample calling, and force call the alleles NOT USING POP PRIORS

    // Scatter sites vcf into multiple chunks
    SCATTER_VCF (
        ch_sites_to_genotype,
        params.variants_per_chunk
    )

    SCATTER_VCF.out.interval_vcf
        .flatMap { vcfs, tbis  ->
            // normalize to a list for cases where there are only 1 bed output for a sample
            def vcfList = (vcfs instanceof List) ? vcfs : [vcfs]
            def tbiList = (tbis instanceof List) ? tbis : [tbis]

            // emit one tuple per bed file
            (0..<vcfList.size()).collect { i ->
                def vcf = vcfList[i] as Path
                def tbiPath = tbiList[i]
                def base = vcf.getFileName().toString()
                base = base.replaceFirst(/\.gz$/, '')
                base = base.replaceFirst(/\.vcf$/, '')
                def interval_hash = base.startsWith('_') ? base.substring(1) : base
                tuple(interval_hash, vcf, tbiPath)
            }
        }
        .filter { interval_hash, vcf, tbi -> vcf && tbi.size() > 0 }   // drop empty
        .set { ch_scatter_vcf }


    // JOINT CALLING WHOLE COHORT 

    // combine sample-level cram with each interval_bed file and interval chunk
    //ch_sample_cram 
    //    .combine ( ch_scatter_vcf )
    //    .map { sample, cram, crai, interval_hash, vcf, tbi -> [ interval_hash, cram, crai ] }
    //    .groupTuple ( by: [0] )
    //    // join to get back interval_file
    //    .join ( ch_scatter_vcf, by: [0] )
    //    // variant type and interval hash columns are combined into a single string for compatibility with mpileup
    //    .map { interval_hash, cram, crai, vcf, tbi -> tuple(interval_hash, vcf, tbi, cram, crai) }
	//   .set { ch_cram_to_genotype }

    // Calculate cohort size for memory scaling
    //ch_cohort_size = ch_sample_names.unique().count()

    // Call just target sites using mpileup
    //MPILEUP (
    //    ch_cram_to_genotype,
    //    ch_genome_indexed,
    //    ch_cohort_size
    //)
     
    // SINGLE SAMPLE CALLING
    ch_sample_cram 
        .combine ( ch_scatter_vcf )
        .map { sample, cram, crai, interval_hash, vcf, tbi -> tuple(sample, interval_hash, vcf, tbi, cram, crai) }
        .set { ch_cram_to_genotype }

    // Call just target sites using mpileup
    MPILEUP_SINGLE (
        ch_cram_to_genotype,
        ch_genome_indexed
    )

    // Merge seperate VCFs by sample
    MPILEUP_SINGLE.out.vcf
        .groupTuple(by: 0)
        .set { ch_vcf_to_merge }

    MERGE_VCFS (
        ch_vcf_to_merge
    )

    emit: 
    vcf = MPILEUP.out.vcf

}