process PLINK_PCA {
    def process_name = "plink_pca"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/plink", mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "PLINK/2.00a3.7-gfbf-2024a"

    input:
    tuple val(outname), path(plinkfiles)

    output: 
    tuple val(outname), path("*.pca*"),                           emit: pca
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash

    # Export PLINK2 parameters
    export MIN_MAC='${params.min_mac}'
    export MIN_MAF='${params.min_maf}'
    export SITE_MAX_MISSING='${params.site_max_missing}'
    export SAMPLE_MAX_MISSING='${params.sample_max_missing}'
    export LD_WINSIZE='${params.ld_winsize}'
    export LD_STEPSIZE='${params.ld_stepsize}'
    export LD_THRESHOLD='${params.ld_threshold}'
    export PLINK_NUM_PC='${params.plink_num_pc}'

    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.mega} \
        "${outname}"
    """
}