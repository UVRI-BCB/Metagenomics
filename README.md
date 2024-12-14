# Metagenomics analysis

Basic steps of the analysis pipeline

1. Quality check of raw data - Fastqc/MultiQC
2. Reads are quality trimmed and filtered using trim-galore
3. Taxonomic identified of trimmed reads using Kraken2
4. Human reads are removed using Bowtie2.
5. SPADES/IDBA - are used to perform de novo assembly on the human-free reads.
6. The resulting contigs are mapped to a set of specific reference genomes using minimap2.
7. The best matching reference for each contig is selected for short read mapping.
8. The human-free reads from step 4 are aligned to the selected reference genomes using minimap2
9. Sequence variants are called from the alignments using bcftools
10. Sequence variants are applied to the corresponding reference sequences to create consensus sequences.
