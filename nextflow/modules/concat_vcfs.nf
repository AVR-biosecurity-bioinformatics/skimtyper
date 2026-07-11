process CONCAT_VCFS {
    def process_name = "concat_vcfs"    
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    module "BCFtools/1.21-GCC-13.3.0"

    input:
    tuple val(outname), path(vcf), path(vcf_tbi)
    
    output: 
    tuple val(outname),  path("${outname}.{vcf,g.vcf}.gz"), path("${outname}.{vcf,g.vcf}.gz.tbi"),       emit: vcf
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
     
    # Write list of vcf files to process
    printf "%s\n" ${vcf} | sort > vcf.list

    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        "${outname}"

    """
}