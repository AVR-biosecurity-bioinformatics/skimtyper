process SCATTER_VCF {
    def process_name = "scatter_vcf"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BCFtools/1.22-GCC-13.3.0"

    input:
    tuple path(vcf), path(tbi)
    val(counts_per_chunk)

    output: 
    tuple path("*.vcf.gz"), path("*.vcf.gz.tbi"),    emit: interval_vcf
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
   
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        ${vcf} \
        ${counts_per_chunk} 

    """
  
}