

##########################################################################################
# Libraries
library(dndscv)
library(Rsamtools)
library(GenomicRanges)
library(Biostrings)

complement = vector()
complement["A"] = "T";	complement["C"] = "G"
complement["G"] = "C";	complement["T"] = "A"

args = commandArgs(trailingOnly=TRUE)
out_prefix    = args[1]
mut_file      = args[2]
tbix_file     = args[3]
regions_file  = args[4]

base_dir  = "."

muts = read.table(mut_file,sep="\t",header=T)
muts$ref  = sapply(muts$context,function(x) unlist(strsplit(x,""))[2])
muts$pyr1 = sapply(muts$pyrsub, function(x) unlist(strsplit(x,""))[2])
muts$to1  = sapply(muts$pyrsub, function(x) unlist(strsplit(x,""))[5])
muts$pyr_sub = paste(sapply(muts$pyrsub, function(x) unlist(strsplit(x,""))[2]),">",sapply(muts$pyrsub, function(x) unlist(strsplit(x,""))[5]),sep="")
complement_these = which(muts$ref != muts$pyr1)
muts$to = muts$to1
muts[complement_these,"to"] = complement[muts[complement_these,"to1"]]
muts$sub = paste(muts$ref,">",muts$to,sep="")

regions           = read.table(paste(base_dir,"/",regions_file,sep=""),sep="\t",header=F)[,1:3]
colnames(regions) = c("chr","start","end")
regions$start     = regions$start+1 # 1-based
regions$track     = regions_file
regions$feature   = "na"
regions$chr       = gsub("chr","",regions$chr)
regions           = regions[which(regions$chr %in% c(c(1:22),"X"),]


##########################################################################################
# Get sequence composition of each region:
# Use TabixFile from Rsamtools package
tbx=TabixFile(tbix_file)

sub_vec        = c("C>A","C>G","C>T","T>A","T>C","T>G")
ctx_vec        = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
full_vec       = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
xstr           = paste(substr(full_vec,5,5), substr(full_vec,1,1), substr(full_vec,7,7), sep="")
ordered_names  = paste(xstr,">",rep(c("A","G","T","A","C","G"),each=16),sep="")
tris           = unique(xstr)
tris_compl     = vector()
for(tri in tris) {
	tt = unlist(strsplit(tri,""))
	tris_compl[tri] = paste(complement[tt[3]],complement[tt[2]],complement[tt[1]],sep="")
}

get_base_counts = function(tbx,regions) {
	res <- scanTabix(tbx, param=regions)
	lengths = sapply(res,length)
	counts = matrix(nrow=length(regions),ncol=length(tris)*2)
	colnames(counts) = c(tris,tris_compl)
	for(i in 1:length(regions)) {
		if(lengths[i]==0) {
			counts[i,] = 0
		} else {
			res2 = t(sapply(res[[i]],function(x) unlist(strsplit(x,"(\t)|(;)"))[c(4,6)]))			
			counts[i,] = sapply(colnames(counts),function(x) sum(as.numeric(res2[which(res2[,1]==x),2])))
		}
	} 
	counts[,tris] = counts[,tris]+counts[,tris_compl]
	counts = counts[,tris]
	return(counts)
}

all_counts = as.data.frame(matrix(nrow=0,ncol=5+length(tris)))
colnames(all_counts) = c("track","feature","chr","start","end",tris)

start = 1
size  = 500
while(start < nrow(regions)) {
	cat(start,"/",nrow(regions),"\n")
	end = min(nrow(regions),start+size-1)
	param = GRanges(regions[start:end,"chr"], IRanges(start=regions[start:end,"start"]-1,end=regions[start:end,"end"]))
	counts = get_base_counts(tbx,param)
	j = regions[start:end,c("track","feature","chr","start","end")]
	j[,colnames(counts)] = counts
	if(nrow(all_counts) == 0) {
		all_counts = j
	} else {
		all_counts = rbind(all_counts,j)
	}
	start = start+size
}


# Intersect mutations:
muts_gr        <- GRanges(muts$chrom,         IRanges(start=muts$chromStart,         end=muts$chromStart  ))
regions_gr     <- GRanges(all_counts[,"chr"], IRanges(start=all_counts[,"start"]-1, end=all_counts[,"end"]))
muts_in_regions = as.data.frame(findOverlaps(muts_gr,regions_gr))

for(pyrsub in ordered_names) {
	all_counts[,pyrsub] = 0
}

if(nrow(muts_in_regions) > 0) {
	for(i in 1:nrow(muts_in_regions)) {
		sub = muts[muts_in_regions$queryHits[i],"pyrsub"]
		all_counts[muts_in_regions$subjectHits[i],sub] = all_counts[muts_in_regions$subjectHits[i],sub] + 1
	}
}
all_counts = apply(all_counts[,c(tris,ordered_names)],2,sum)
kk=c(regions_file,epi_feat,out_prefix,all_counts)
names_ = names(kk)
names_[1:3] = c("regions_file","epi_feat","out_prefix")
names(kk) = names_
write.table(t(kk),paste(out_prefix,".",regions_file,".trinuc.summary.tsv",sep=""),sep="\t",row.names=F,col.names=T,quote=F)


