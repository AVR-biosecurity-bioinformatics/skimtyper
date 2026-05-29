process PLINK_IMPORT {
    def process_name = "plink_import"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/plink", mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "PLINK/2.00a3.7-gfbf-2024a"

    input:
    tuple val(outname), path(vcf), path(vcf_tbi)

    output: 
    path("${outname}.{bim,bed,fam}"),                           emit: plink
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash

    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        "${outname}" \
        ${vcf}
    """
}