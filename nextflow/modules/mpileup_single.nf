process MPILEUP_SINGLE {
    def process_name = "mpileup"    
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BCFtools/1.21-GCC-13.3.0:BEDTools/2.31.1-GCC-13.3.0:SAMtools/1.22.1-GCC-13.3.0:parallel/20240722-GCCcore-13.3.0"

    input:
    tuple val(sample), val(interval_hash), path(vcf), path(tbi), path(cram), path(cram_index)
    tuple path(ref_genome), path(genome_index_files)

    output: 
    tuple val(sample), val(interval_hash), path("${sample}.${interval_hash}.vcf.gz"), path("${sample}.${interval_hash}.vcf.gz.tbi"),    emit: vcf

    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    # Export Mpileup parameters
    export RMDUP='${params.rmdup}'
    export EXCLUDE_PAD='${params.exclude_padding}'
    export OUTPUT_INVARIANT='${params.output_invariant}'
    export PLOIDY='${params.ploidy}'
    export MINBQ='${params.minbq}'
    export MINMQ='${params.minmq}'
    export MIN_ALIGNED_LENGTH='${params.min_aligned_length}'
    export MIN_FRAGMENT_LENGTH='${params.min_fragment_length}'
    export MAX_FRAGMENT_LENGTH='${params.max_fragment_length}'
    export MUTATION_RATE='${params.mutation_rate}'
    export MAXDEPTH='${params.max_depth}'

    # Write list of cram files to process
    printf "%s\n" ${cram} | LC_ALL=C sort -u > cram.list

    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        ${ref_genome} \
        ${sample} \
        ${interval_hash} \
        ${vcf} 

    """
}