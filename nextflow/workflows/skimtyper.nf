

//// import subworkflows
include { VALIDATE_INPUTS                                           } from '../subworkflows/validate_inputs'
include { PROCESS_READS                                             } from '../subworkflows/process_reads'
include { GENOTYPE_WITH_PANEL                                       } from '../subworkflows/genotype_with_panel'

//// import modules
include { INDEX_GENOME                                              } from '../modules/index_genome' 
include { INDEX_MITO                                                } from '../modules/index_mito'


// Create default channels
ch_dummy_file = file("$baseDir/assets/dummy_file.txt", checkIfExists: true)
ch_reports = Channel.empty()
ch_multiqc_config   = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)

workflow SKIMTYPER {

    /*
    Input channel parsing
    */    

    if ( params.samplesheet ){
        ch_samplesheet = Channel
            .fromPath (
                params.samplesheet,
                checkIfExists: true
            )
    } else {
        println "\n*** ERROR: 'params.samplesheet' must be given ***\n"
    }
    
    // Parse input samplesheet
    ch_samplesheet
        .splitCsv(header: true)
        .map { row ->
            // Fail early if required columns are missing
            def required = ['sample','pop','fwd','rev']
            def present  = row.keySet()*.toString() as Set
            def missing  = required.findAll { !(it in present) }
            if( missing ) {
                error "Samplesheet is missing required columns: ${missing.join(', ')}. " +
                    "Found columns: ${present.toList().sort().join(', ')}"
            }

            // Parse samplesheet columns
            def sample = row.sample.toString().trim()
            def pop    = row.pop.toString().trim().replaceAll(/\s+/, '_')
            def r1     = file(row.fwd, checkIfExists: true)
            def r2     = file(row.rev, checkIfExists: true)

            // Fail early if any sample names less than 3 characters
            if( sample.size() < 3 ) {
                error "Invalid sample name '${sample}' in samplesheet. " +
                    "Sample names must be at least 3 characters long because bcftools +fill-tags fails on 2-character sample IDs."
            }

            def lib = r1.getName().replaceFirst(/\.(fastq|fq)\.gz$/, '')
            tuple(sample, lib, pop, r1, r2)
        }
        .set { ch_samplesheet_parsed }

    // Reads channel
    ch_samplesheet_parsed
        .map { sample, lib, pop, r1, r2 -> tuple(sample, lib, r1, r2) }
        .set { ch_reads }

    // Sample names channel
    ch_samplesheet_parsed
        .map { sample, lib, pop, r1, r2 -> sample }
        .unique()
        .set { ch_sample_names }

    // Sample names and pops channel
    ch_samplesheet_parsed
        .map { sample, lib, pop, r1, r2 -> tuple(sample, pop) }
        .set { ch_sample_pop }

    // Reference genome channel
    if ( params.ref_genome ){
        ch_genome = Channel
            .fromPath (
                params.ref_genome, 
                checkIfExists: true
            )
    } else {
        ch_genome = Channel.empty()
    } 
    

    // Reference panel channel: VCF + tabix index
    ch_panel_vcf = Channel.fromPath(
        params.ref_panel,
        checkIfExists: true
    )

    ch_panel_tbi = Channel.fromPath(
        "${params.ref_panel}.tbi",
        checkIfExists: true
    )

    ch_panel_vcf
        .combine(ch_panel_tbi)
        .map { vcf, tbi -> tuple("panel", vcf, tbi) }
        .set{ ch_panel}

    /*
    Process nuclear genome
    */

    INDEX_GENOME (
        ch_genome
    )

    ch_genome_indexed = INDEX_GENOME.out.fasta_indexed.first()
    ch_genome_bed = INDEX_GENOME.out.genome_bed
    
    
    /*
    Process mitochondrial genome and create intervals
    */
        
    INDEX_MITO (
        ch_genome,
        params.mito_contig
    )

    ch_mito_indexed = INDEX_MITO.out.fasta_indexed.first()
    ch_mito_bed = INDEX_MITO.out.bed.first()
    
    /*
    Validate inputs
    */

    VALIDATE_INPUTS (
        ch_sample_names,
        ch_reads,
        ch_genome_indexed
    )

    /*
    Process reads per sample, aligning to the genome, and merging
    */

    PROCESS_READS (
        ch_sample_names,
        VALIDATE_INPUTS.out.validated_fastq,
        VALIDATE_INPUTS.out.rg_to_validate,
        ch_genome_indexed
    )
    
    PROCESS_READS.out.perbase
        .set{ ch_read_counts }

    /*
    Process reads per sample, aligning to the genome, and merging
    */
     GENOTYPE_WITH_PANEL (
        ch_sample_names,
        PROCESS_READS.out.cram,
        ch_genome_indexed,
        ch_panel,
        ch_read_counts
    )

    GENOTYPE_WITH_PANEL.out.vcf
        .set{ ch_vcfs }

}