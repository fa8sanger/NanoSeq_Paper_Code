# NanoSeq_Paper_Code

This repository contains the code for the main downstream analyses presented in the 
manuscript. The repository for  processing of NanoSeq libraries, calling substitutions, and calculating
mutation burden is available at:

https://github.com/cancerit/NanoSeq


Below is a brief summary of the content of each directory.

### DATA
Data is in the DATA directory. 

### INDEL_CALLING
Collection of Perl and R scripts to call indels in NanoSeq data. Instructions are provided in indelCalling_pipeline.sh

### FIGURE_1, FIGURE_2, FIGURE_3 
Contain the code to reproduce the main figures, and carry the analyses involved. 

Script sections refer to specific figure panels. 

Code for panels 3f and 3g has not been included because dependencies with third-party data files and multiple, heavy data files.  

The scripts also contain code for related analyses. For instance, FIGURE_3/Figure3.R contains the code involved in de novo signature extraction; FIGURE_2/Figure2.R contains the code to do the lineal regression of blood data.

### BARCODE_CLASHES
Contains R code to estimate the impact of barcode clashes (chimeric read bundles) on burden estimates

### BURDEN_IN_SPECIFIC_REGIONS
Code to estimate burden in specific genomic regions such as all genes, 
expressed genes or heterochromatin regions. Instructions are provided in the shell
script named burden_specific_regions_pipeline.sh 

### CORRECTION_BURDEN_STD_SEQUENCING
Code and data to estimate mutation burden in standard sequencing data.
Files human_genome.bed and nanoseq_genome.bed are available at https://drive.google.com/drive/folders/1wqkgpRTuf4EUhqCGSLA4fIg9qEEw3ZcL?usp=sharing

### FMOL_VS_LIBCOMPLEXITY
Code to estimate the relationship between library yields (in fmol) and library complexity (in number of molecules).

### GENOME_MASKS
Genomic masks for common SNP masking and detection of noisy/variable genomic sites. Build GRCh37.
These masks are available at: https://drive.google.com/drive/folders/1wqkgpRTuf4EUhqCGSLA4fIg9qEEw3ZcL?usp=sharing

### NEURONS_SINGLECELLSEQ_SIGNATURES
R code and data for analysing Lodato et al (Nature 2018) single cell sequencing data, including signature extraction and signature comparison.

### SIGNATURE_RATES
This code allows studying accumulation of mutations with age for each signature, both for neurons and smooth muscle. 
The code is associated with Exteded Data Figure 7, panels e and f.


