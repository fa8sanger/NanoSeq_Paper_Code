##########################################################################################
# Script to correct burden estimates based on minimum coverage thresholds (calculated with
#  min_cov_genome.pl) and coverage in regions passing NanoSeq filters
# Trinucleotide frequencies correction are also applied
# Modify the genome_file variable accordingly
#
# As currrently implemented it expects a Caveman output file with name:
#  [sample].caveman_c.flag.vcf.gz
# Modify to your requirements. This file is only needed to get the list of mutations and
# their coordinates.
# The script also expects to find a file named [sample].min20.bed generated with
#  the script "min_cov_genome.pl"
##########################################################################################

library(GenomicRanges)
library(deepSNV)
setwd("CORRECTION_BURDEN_STD_SEQUENCING")
genomeFile             = "/lustre/scratch116/casm/cgp/users/fa8/data/hs37d5.fa"

args                   = commandArgs(TRUE) 
c20_genome             = args[1]
nano_genome            = "nanoseq_genome.bed"
     
tri_corrections_       = read.table("trint_counts_and_ratio2genome.tsv",sep="\t",header=T,row.names=1)
tri_corrections        = tri_corrections_[,2]
names(tri_corrections) = rownames(tri_corrections_)

sample = c20_genome

##########################################################################################
# Find the intersection between the C20 genome and the NanoSeq genome:
c20_genome = paste(sample,".min20.bed",sep="")
cat("Loading nanoseq genome...\n")
nano = read.table(nano_genome,sep="\t",header=F)
cat("Loading c20 genome...\n")
c20  = read.table(c20_genome,sep="\t",header=F)
cat("Intersections...\n")
nano_gr = GRanges(nano[,1], IRanges(nano[,2]+1, nano[,3]))
c20_gr  = GRanges(c20[,1],  IRanges(c20[,2]+1,  c20[,3]))

##########################################################################################
# Load mutations:
cat("Loading muts...\n")
muts_file = paste(sample,".caveman_c.flag.vcf.gz",sep="")
zz=gzfile(muts_file,'rt')   
muts_caveman=read.table(zz,header=F)
muts = muts_caveman[which(muts_caveman[,7] == "PASS"),]

##########################################################################################
# Find overlaps between mutations and regions:
muts_gr = GRanges(muts[,1], IRanges(muts[,2], muts[,2]))
nano_c20_gr = GenomicRanges::intersect(nano_gr, c20_gr)
nano_c20 = as.data.frame(nano_c20_gr)

sum(as.numeric(c20[,3]-c20[,2]))
sum(as.numeric(nano[,3]-nano[,2]))
sum(as.numeric(nano_c20[,3]-nano_c20[,2]))

muts_in_nano_c20 = as.matrix(findOverlaps(muts_gr,nano_c20_gr))
muts_in_c20      = as.matrix(findOverlaps(muts_gr,c20_gr))

mut_indexes = muts_in_nano_c20[,1]
final_muts  = muts[mut_indexes,]

##########################################################################################
# Get trinucleotide/pyrimidine context for each mutation. Necessary to apply the trint
# composition corrections
trints = as.vector(scanFa(genomeFile, GRanges(final_muts[,1], 
                                      IRanges(final_muts[,2]-1, final_muts[,2]+1))))
complement = vector()
complement["A"] = "T"; complement["C"] = "G"; complement["G"] = "C"; complement["T"] = "A";		
for(i in c(1:length(trints))) {
	trint = trints[i]
	tri = unlist(strsplit(trint,""))
	if(tri[2] %in% c("G","A")) { # reverse complement if ref base is a purine
		trint = paste(complement[tri[3]],complement[tri[2]],complement[tri[1]],sep="")
		trints[i] = trint
	}
}

##########################################################################################
# Output corrected mutation rates:
cat("RESULT\t",muts_file,
          "\t",sample,
          "\t",nrow(muts),                                                                    # original number of mutations
          "\t",2*sum(as.numeric(c20[,3]-c20[,2])),                                            # coverage in the C20 genome (diploid)
          "\t",nrow(muts_in_c20),"\t",nrow(muts_in_nano_c20),                                 # number of mutations in C20
          "\t",sum(1/tri_corrections[trints]),                                                # number of mutations in C20+NanoSeq corrected by trinucleotide composition
          "\t",2*sum(as.numeric(nano_c20[,3]-nano_c20[,2])),                                  # coverage in the C20+NanoSeq genome (diploid)
          "\t",nrow(muts_in_nano_c20)/(2*sum(as.numeric(nano_c20[,3]-nano_c20[,2]))),         # burden in the C20+NanoSeq genome
          "\t",sum(1/tri_corrections[trints])/(2*sum(as.numeric(nano_c20[,3]-nano_c20[,2]))), # burden in the C20+NanoSeq genome corrected by trinucleotide composition
          "\n", sep="")


##########################################################################################

