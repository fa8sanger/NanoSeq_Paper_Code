library(deepSNV)
library(vcfR)
library("GenomicRanges")
library("Rsamtools")
library("MASS")

FLANK = 5
genomeFile = "/lustre/scratch116/casm/cgp/users/fa8/data/hs37d5.fa"



args = commandArgs(trailingOnly=TRUE)
vcf_file = args[1]
bam_file = args[2]

vcf <- read.vcfR( vcf_file, verbose = FALSE )
vcf_mat <- data.frame(CHROM=getCHROM(vcf),POS=getPOS(vcf),ID=getID(vcf),REF=getREF(vcf),ALT=getALT(vcf),QUAL=getQUAL(vcf),FILTER=getFILTER(vcf),INFO=getINFO(vcf))
vcf = vcf_mat

# Create regions:
#  * For deletions: get pos + length(deletion)   // length(deletion) = length(ref)-1
#  * For insertions: get pos
#  * Then, sum +/-5 to each side to count indels in the vecinity


for(i in c(1:nrow(vcf))) {
	pos = vcf[i,"POS"]
	chr = vcf[i,"CHROM"]
	len = max(1,length(vcf[i,"REF"])-1)
	start = pos - FLANK
	end   = pos + len + FLANK
	kk = bam2R(bam_file, chr,start,end, q=-100, mask=3844, mq=10) 
	# dont want a filter on BQ because in some bams BQ of indels have -1
    n_bases  = sum(kk[,c("A","C","G","T","a","c","g","t")])
	n_indels = sum(kk[,c("-","INS","_","ins")            ]) # Number of reads with an indel around the mutation
	cat(chr,pos,len,n_bases,n_indels,"\n");
	cat(n_bases,"/",n_indels,"\n")
	if(n_bases == 0) {
		vcf[i,"FILTER"] = paste("MISSINGBULK[",n_indels,"/",n_bases,"]",sep="")
	
	} else if(n_indels/n_bases > 0.01) {
		vcf[i,"FILTER"] = paste("NEI_IND[",n_indels,"/",n_bases,"]",sep="")
	}	
	sequence = as.vector(scanFa(genomeFile, GRanges(chr, IRanges(start-3, end+3))))
	vcf[i,"FILTER"] = paste(vcf[i,"FILTER"] ,";",sequence,sep="")
}

vcf_file = gsub(".vcf",".bulkmasked.vcf",vcf_file)
vcf_file = gsub(".vcf",".tsv",vcf_file)
write.table(vcf,file=vcf_file,row.names=F,col.names=T,sep="\t",quote=F)









