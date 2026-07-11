process PUBLISH_VCF {
    def process_name = "publish_vcf"    
    publishDir "${launchDir}/output/results/vcf", mode: 'copy'

    input:
    tuple val(outname), path(vcf), path(vcf_tbi)
    
    output: 
    tuple path("final.vcf.gz"), path("final.vcf.gz.tbi"), emit: vcf
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
     
    ln -s ${vcf} final.vcf.gz
    ln -s ${tbi} final.vcf.gz.tbi

    """
}
