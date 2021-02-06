# NanoSeq_Paper_Code

This repository contains the code for the main downstream analyses presented in the 
manuscript. The repository for upstream processing NanoSeq libraries, calling substitutions, and calculating
mutation burden is available at:

https://github.com/cancerit/NanoSeq


Below is a brief summary of the content of each directory.

### DATA
Contains... data. 

### INDEL_CALLING
Collection of Perl and R scripts to call indels in NanoSeq data. Instructions are provided in indelCalling_pipeline.sh
Eventually this may be incorporated into the upstream processing pipeline (https://github.com/cancerit/NanoSeq)

### FIGURE_1, FIGURE_2, FIGURE_3 
Contain the code to reproduce the three main figures in the manuscript and related analyses. 

Script sections refer to specific figure panels. 

The scripts also contain code for associated analyses and results presented in the text. For instance, FIGURE_3/Figure3.R contains the code involved in de novo signature extraction; or FIGURE_2/Figure2.R contains the code to do different lineal regressions of blood data.

Code for panels 3f and 3g has not been included because of dependencies with third-party data files and multiple, heavy, downstream data files. Links to third-party data are provided in FIGURE_3/Figure3.R. Instructions to generate the required NanoSeq data files are provided in BURDEN_IN_SPECIFIC_REGIONS.  

### BURDEN_IN_SPECIFIC_REGIONS
Code to estimate burden in specific genomic regions such as in each gene, in groups of genes or in heterochromatin regions. Instructions are provided in the shell script named burden_specific_regions_pipeline.sh 

### BARCODE_CLASHES
Contains R code to estimate the impact of barcode clashes (chimeric read bundles) on burden estimates

### FMOL_VS_LIBCOMPLEXITY
Code to estimate the relationship between library yields (in fmol) and library complexity (in number of molecules).

### CORRECTION_BURDEN_STD_SEQUENCING
Code and data to estimate mutation burden in standard sequencing data.
Files human_genome.bed and nanoseq_genome.bed are available at https://drive.google.com/drive/folders/1wqkgpRTuf4EUhqCGSLA4fIg9qEEw3ZcL?usp=sharing

### GENOME_MASKS
Genomic masks for common SNP masking and detection of noisy/variable genomic sites. Build GRCh37.
These masks are available at: https://drive.google.com/drive/folders/1wqkgpRTuf4EUhqCGSLA4fIg9qEEw3ZcL?usp=sharing

### NEURONS_SINGLECELLSEQ_SIGNATURES
R code and data to analyse Lodato et al (Nature 2018) single cell sequencing data, including signature extraction and signature comparison.

### SIGNATURE_RATES
This code allows studying accumulation of mutations with age for each signature, both for neurons and smooth muscle. 
The code is associated with Exteded Data Figure 7, panels e and f.


