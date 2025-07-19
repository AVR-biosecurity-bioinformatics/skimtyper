# `clueseq`
A nextflow-based pipeline for rapid and genotype free relatedness and sample origin prediction from genomic data.

> This pipeline is currently **Experimental** and being actively developed, with no guarantee that the code is stable or usable!

`clueseq` aims to rapidly predict relatedness between query samples and predict their geographic and population origin, in order to provide 'clues' for time-sensitive genomic investigations of biosecurity outbreaks.

`clueseq` requires as inputs:
- A set of query samples in fastq format
- A set of reference samples in fastq or BAM format
- A metadata sheet annotating the population and geographic coordinates of the reference samples
- A panel of known reliable variant positions used for extracting 'genotype-like' information

`clueseq` is a sibling pipeline to [`skimseq`](https://github.com/AVR-biosecurity-bioinformatics/skimseq), which conducts traditional traditional variant calling and can be used to generate a list of reliable variant positions for use with `clueseq`.

Pipeline ideas:
- Input read QC and filtering with `fastp`
- Reference alignment with `BWA-mem2`
- Alternatively, alignment to just informative positions
- Extraction of 'genotype-like' data using [`somalier`](https://github.com/brentp/somalier)
- Relatedness estimation using [`somalier`](https://github.com/brentp/somalier)
- Supervised PCA using [`somalier`](https://github.com/brentp/somalier)
- Geographic placement using [`ReLocator`](https://github.com/kr-colab/ReLocator)
