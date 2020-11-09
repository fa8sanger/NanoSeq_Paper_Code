# Federico Abascal, 2020
# Pipeline (bash script) to call indels from NanoSeq data 
# Requirements: Perl, R, samtools, bcftools
# Required R libraries: deepSNV, vcfR, GenomicRanges, Rsamtools and MASS

# 1a. Change to the directory containing the tables produced by the NanoSeq pipeline
cd [sampleX_directory_tables]

# 1b. Create a directory to write temporary indel calling files
mkdir ../sampleX_indels


# 1c. Run first part of the pipeline to preprocess the tables and identify read bundles of interest:
#     This can be paralelised for each individual table in [sampleX_directory_tables]
#     Parameters indicate: 
#       1. minimum number of reads per strand in read bundle, 
#       2. number of bases to trim from the 5' end of reads,
#       3. and number of bases to trim from 3' end of reads
for a in *bed.gz; 
   do zcat $a | perl indelCaller_step1.pl 2 10 135 > ../sampleX_indels/$a.indels.bed
done


# 2a. Change to the directory containing the temporary indel calling files:
cd ../sampleX_indels

# 2b. Run the second step of the pipeline, which is the main indel caller, and relies on 
#     samtools and bcftools. 
#     Two arguments are required:
#       1. The BAM file containing the NanoSeq sequence data
#       2. A prefix for the output file name
#     Modify the path to the human reference genome (line #19)
#     This perl script requires bcftools (tested with version 1.9) and samtools 
#       (tested with version 1.9) are in the path.
#     Indels overlapping common SNPs or noisy sites (see our SNP and noise masks) are
#       flagged as 'MASKED'
grep -hv BULK *bed | perl indelCaller_step2.pl sampleX_nanoseq.bam sampleX.indels


# 3. Finally, run the following Rscript to check whether each of the potential indel calls
#    happen at loci rich in indels in the matched normal. Those would be flagged as 'NEI_IND'.
#    Reliable indels are flagged with 'PASS'.
#    Two parameters are expected:
#      1. The name of the output file from step #5 [sampleX.indels.final.vcf]
#      2. The BAM file containing the matched normal
#    Modify line #8 to change the path to the human reference genome
Rscript indelCaller_step3.R sampleX.indels.final.vcf matched_normal.bam


