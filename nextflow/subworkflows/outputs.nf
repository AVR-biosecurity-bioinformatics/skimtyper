/*
    Create outputs
*/

//// import modules
include { VCF2DIST                                               } from '../modules/vcf2dist' 
include { PLOT_ORDINATION                                        } from '../modules/plot_ordination' 
include { PLOT_TREE                                              } from '../modules/plot_tree' 
include { PLINK_IMPORT                                           } from '../modules/plink_import' 
include { PLINK_PCA                                              } from '../modules/plink_pca' 

workflow OUTPUTS {

    take:
    ch_vcf
    ch_sample_pop

    main: 

    /* 
        Create outputs
    */

    // Create distance matrices from VCFs
    VCF2DIST (
        ch_vcf
    )


    // Import PLINK file
    PLINK_IMPORT (
        ch_vcf
    )

    // Run PLINK PCA
    PLINK_PCA (
        PLINK_IMPORT.out.plink
    )   

    // Turn ch_sample_pop tuples into a 2‑col TSV 'popmap' file
    ch_sample_pop
        .map { s,p -> "$s\t$p" }
        .collectFile(name: 'sample_pop.tsv', newLine: true)
        .first()
        .set { ch_popmap }

    // create ordination plot from distance matrices
    PLOT_ORDINATION (
        VCF2DIST.out.mat,
        ch_popmap,
        false
    )

    // Create NJ tree
    PLOT_TREE (
        VCF2DIST.out.mat,
        ch_popmap
    )

}