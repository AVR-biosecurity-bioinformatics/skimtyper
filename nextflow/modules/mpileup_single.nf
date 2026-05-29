process MPILEUP_SINGLE {
    def process_name = "mpileup_single"    
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BCFtools/1.23.1-GCC-13.3.0"

    input:
    tuple val(sample), val(interval_hash), path(panel_vcf), path(panel_tbi), path(cram), path(cram_index)
    tuple path(ref_genome), path(genome_index_files)

    output: 
    tuple val(sample), val(interval_hash), path("${sample}.${interval_hash}.vcf.gz"), path("${sample}.${interval_hash}.vcf.gz.tbi"),    emit: vcf

    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    # Export Mpileup parameters
    export RMDUP='${params.rmdup}'
    export PLOIDY='${params.ploidy}'
    export MINBQ='${params.minbq}'
    export MINMQ='${params.minmq}'
    export MIN_ALIGNED_LENGTH='${params.min_aligned_length}'
    export MIN_FRAGMENT_LENGTH='${params.min_fragment_length}'
    export MAX_FRAGMENT_LENGTH='${params.max_fragment_length}'
    export MUTATION_RATE='${params.mutation_rate}'
    export MAXDEPTH='${params.max_depth}'
    export REF_PANEL_PRIORS='${params.ref_panel_priors}'

    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        ${ref_genome} \
        ${sample} \
        ${interval_hash} \
        ${panel_vcf} \
        ${cram}

    """
}