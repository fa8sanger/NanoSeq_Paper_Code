# Federico Abascal, 2020
# Pipeline (bash script) to estimate the mutation rate in specific regions using NanoSeq data 
# This pipeline takes into account the uneven coverage of NanoSeq data and the variable
#   trinucleotide frequencies.
# Requirements: NanoSeq pipeline, gzip, tabix, bedtools, awk
# Required R libraries: dndscv, Rsamtools, GenomicRanges, Biostrings



# 1a. To run this pipeline, the default NanoSeq pipeline needs to be modified.
#     When running the variantcaller program, option -U has to be invoked, providing a file
#     name where site-specific coverage is written.
#     Because these files can be very large, gzip compression is recommended.
#     The code below is to be run in parallel for each specific output table.
variantcaller  [...] [ other parameters ] -U outdir/sampleX.coverage; 
gzip -f $outdir/sampleX.coverage

# 1b. Once the coverage files for all tables have been generated, we will merge them into
#     a single, compressed and indexed file
zcat variantcalls-sampleX/*coverage.gz | uniq -c | awk '{print \$2,\$3,\$4,\$5\";\"\$1;}' | \
     tr ' ' '\t' | sort -k 1,1 -k2,2n -T . -S 2G | bgzip > variantcalls-$tag/genome.coverage.bed.gz; \
     tabix -p bed variantcalls-$tag/genome.coverage.bed.gz
     
# 2. Run Rscript to calculate coverage in the regions of interest
#    Arguments:
#      1. Output prefix name
#      2. Mutation file, this is the output file of the NanoSeq pipeline which ends in "muts.tsv"
#      3. Tabix file, this is the genome.coverage.bed.gz file
#      4. Regions file, this is a bed file containing the coordinates of the regions of interest
#         These regions may be highly transcribed genes, heterochromatin regions, etc
#    The results of running this script is a file named [chosen output prefix name].trinuc.summary.tsv
#    This file contains the coverage for each of the 32 pyrimidine trinucleotides and the
#      count of each of 96 trinucleotide substitutions.
Rscript botseq_coverage_and_composition_trinuc_version.R sampleX_in_genes sampleX.muts.tsv \
        variantcalls-sampleX/genome.coverage.bed.gz gene_coordinates.bed

# 3. With the output provide by the previous script, the burden and the substitution profiles 
#    can be obtained with the following R script:
Rscript final_processing.R sampleX_in_genes.trinuc.summary.tsv sampleX_in_genes

