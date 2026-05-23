# `skimtyper`
A low-fat nextflow-based pipeline for genotyping low coverage genome sequencing (skim sequencing) data using a reference panel VCF

> This pipeline is currently **Experimental** and being actively developed, with no guarantee that the code is stable or usable!

`skimtyper` aims to rapidly genotype skim sequencing data using a reference panel VCF, calculate relatedness between query samples, and predict their geographic and population origin for time-sensitive genomic investigations of biosecurity outbreaks.

`skimtyper` requires as inputs:
- A set of query samples in fastq format
- A reference panel in VCF format containg known relaible variant positions

`skimtyper` is a sibling pipeline to [`skimseq`](https://github.com/AVR-biosecurity-bioinformatics/skimseq), which conducts traditional traditional variant calling and can be used to generate a list of reliable variant positions for use with `skimtyper`.

Pipeline steps:
- Read accuracy calculation with [`fraguracy'](https://github.com/brentp/fraguracy)
- Reference alignment with `BWA-mem2`
- Genotype new samples at sites present in the reference panel
- Harmonise and merge new samples with the reference panel
- Calculate relatedness, PCA, and population of-origin prediction using DAPC


Other ideas:
- Extraction of 'genotype-like' data using [`somalier`](https://github.com/brentp/somalier)
- Relatedness estimation using [`somalier`](https://github.com/brentp/somalier)
- Supervised PCA using [`somalier`](https://github.com/brentp/somalier)
- Geographic placement using [`ReLocator`](https://github.com/kr-colab/ReLocator)
