# NanoSeq_Paper_Code

This repository contains the code for the main downstream analyses presented in the 
manuscript.

The repository for  processing of NanoSeq libraries, calling substitutions, and calculating
mutation burden is available at:

https://github.com/cancerit/NanoSeq

Below is a brief summary of the content of each directory.

###BLOOD_REGRESSIONS
Contains R code and accompanying data to do the regression analysis of
granulocytes and stem cell colonies

###BURDEN_IN_SPECIFIC_REGIONS
Code to estimate burden in specific genomic regions such as
expressed genes or heterochromatin, for example. Instructions are provided in the shell
script named burden_specific_regions_pipeline.sh 

###CORRECTION_BURDEN_STD_SEQUENCING
Code and data to estimate mutation burden in standard sequencing data.

###FMOL_VS_LIBCOMPLEXITY
Code to estimate the relationship between library yields (in fmol) and library complexity (in number of molecules).

###GENOME_MASKS
Genomic masks for common SNP masking and detection of noisy/variable genomic sites. Build GRCh37.

###INDEL_CALLING
Collection of Perl and R scripts to call indels in NanoSeq data. Instructions are provided in indelCalling_pipeline.sh

**NEURONS_SINGLECELLSEQ_SIGNATURES**
R code and data for analysing Lodato et al (Nature 2018) single cell sequencing data, including signature extraction and signature comparison.

**SIGNATURE_EXTRACTION**
R code relying on sigfit and data to extract substitution signatures in NanoSeq data.

