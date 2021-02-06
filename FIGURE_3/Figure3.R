##########################################################################################
# Figure 3 panels
##########################################################################################

library(lme4)
library(lmerTest)
library(deconstructSigs)
genome_size = sum(tri.counts.genome)*2

setwd("FIGURE_3")

##########################################################################################
# Accumulation of substitutions and indels with ageg in neurons and smooth muscle
# Panels a, b, j and k

par(mfrow=c(2,2))
for(cell_type in c("NEURONS","SMOOTH_MUSCLE")) {
	# Subs:
	nanoseq                        = read.table("../DATA/rates.tsv",sep="\t",header=T)
	nanoseq$num_muts_per_cell      = nanoseq$burden     * genome_size
	nanoseq$num_muts_per_cell_lci  = nanoseq$burden_lci * genome_size
	nanoseq$num_muts_per_cell_uci  = nanoseq$burden_uci * genome_size
	if(cell_type == "NEURONS") {
		nanoseq = nanoseq[which(nanoseq$cell_type=="neurons"),]	
	} else {
		nanoseq = nanoseq[which(nanoseq$cell_type %in% c("smuscleU","smuscleC")),]	
	}
	nanoseq$donor = substr(nanoseq$sampleID,1,7)
	
	par(mar=c(5,5,3,3))
	par(mgp=c(2, .7, 0))
	plot(NULL,NULL,pch=20,
		 ylab="Number of substitutions per cell",cex=1,
		 main=cell_type,xlab="Age",
		 frame.plot=F,cex.main=.6,cex.lab=.6,cex.axis=.6,
		 col="firebrick",ylim=c(0,max(nanoseq$num_muts_per_cell_uci)),
		 xlim=c(0,100),las=1)
	if(cell_type=="SMOOTH_MUSCLE") {
		mm = lmer(num_muts_per_cell ~ age + (age-1|donor),data=nanoseq)
		m1 = lmer(num_muts_per_cell ~ age + (age-1|donor) + (age-1|cell_type),data=nanoseq)
		anova(mm,m1)
		coef = fixef(mm)
		library(bootpredictlme4)
		ages = seq(from=0,to=100,by=7)
		nsim = 1000
		rr = sapply(ages,function(x) predict(mm, newdata=data.frame(age=x), re.form=NA, se.fit=TRUE, nsim=nsim)$ci.fit[,1])
		polygon(x=c(ages, rev(ages)), y=c(rr[1,], rev(rr[2,])), col="grey90", border=NA)
		segments(nanoseq$age,nanoseq$num_muts_per_cell_lci,nanoseq$age,nanoseq$num_muts_per_cell_uci,col="darkgray")
		points(nanoseq[which(nanoseq$cell_type=="smuscleU"),"age"],nanoseq[which(nanoseq$cell_type=="smuscleU"),"num_muts_per_cell"],pch=20,col="darkolivegreen")
		points(nanoseq[which(nanoseq$cell_type=="smuscleC"),"age"],nanoseq[which(nanoseq$cell_type=="smuscleC"),"num_muts_per_cell"],pch=20,col="darkolivegreen2")
		abline(coef[1], coef[2], col="darkolivegreen",lty=2)		
		legend("topleft",legend=c("Bladder","Colon"),col=c("darkolivegreen","darkolivegreen2"),pch=15,cex=.6,bty = "n")
	} else {
		all_ages <- data.frame(age = c(0:100))
		ci=predict(lm(num_muts_per_cell~age,data=nanoseq),all_ages,interval="confidence", level=0.95)
		polygon(x=c(all_ages$age, rev(all_ages$age)), y=c(ci[,3], rev(ci[,2])), col="grey90", border=NA)
		segments(nanoseq$age,nanoseq$num_muts_per_cell_lci,nanoseq$age,nanoseq$num_muts_per_cell_uci,col="darkgray")
		points(nanoseq[which(nanoseq$disease_state=="Healthy"),"age"],nanoseq[which(nanoseq$disease_state=="Healthy"),"num_muts_per_cell"],pch=20,col="firebrick")
		points(nanoseq[which(nanoseq$disease_state=="AD"),"age"],nanoseq[which(nanoseq$disease_state=="AD"),"num_muts_per_cell"],pch=20,col="firebrick2")
		coef = lm(nanoseq$num_muts_per_cell~nanoseq$age)$coefficients
		abline(coef[1], coef[2], col="firebrick",lty=2)
		legend("topleft",legend=c("Healthy","Alzheimer's disease"),col=c("firebrick","firebrick2"),pch=15,cex=.6,bty = "n")
	}
	
	# Indels:
	nanoseq = nanoseq[which(nanoseq$exclude_indels == "No"),]
	nanoseq$indels_lci      = sapply(nanoseq$indels,function(x) poisson.test(x)$conf.int[1])
	nanoseq$indels_uci      = sapply(nanoseq$indels,function(x) poisson.test(x)$conf.int[2])
	nanoseq$indel_rate      = nanoseq$indels/nanoseq$total
	nanoseq$indel_rate_lci  = nanoseq$indels_lci/nanoseq$total
	nanoseq$indel_rate_uci  = nanoseq$indels_uci/nanoseq$total
	nanoseq$num_indels_per_cell = nanoseq$indel_rate * genome_size
	nanoseq$num_indels_per_cell_lci = nanoseq$indel_rate_lci * genome_size
	nanoseq$num_indels_per_cell_uci = nanoseq$indel_rate_uci * genome_size

	par(mar=c(5,5,3,3))
	par(mgp=c(2, .7, 0))
	plot(NULL,NULL,pch=20,
		 ylab="Number of indels per cell",cex=1,
		 main=cell_type,xlab="Age",
		 frame.plot=F,cex.main=.6,cex.lab=.6,cex.axis=.6,
		 col="firebrick",ylim=c(0,max(nanoseq$num_indels_per_cell_uci)),
		 xlim=c(0,100),las=1)
	if(cell_type=="SMOOTH_MUSCLE") {
		mm = lmer(num_indels_per_cell ~ age + (age-1|donor),data=nanoseq)
		coef = fixef(mm)
		library(bootpredictlme4)
		ages = seq(from=0,to=100,by=7)
		nsim = 1000
		rr = sapply(ages,function(x) predict(mm, newdata=data.frame(age=x), re.form=NA, se.fit=TRUE, nsim=nsim)$ci.fit[,1])
		polygon(x=c(ages, rev(ages)), y=c(rr[1,], rev(rr[2,])), col="grey90", border=NA)
		segments(nanoseq$age,nanoseq$num_indels_per_cell_lci,nanoseq$age,nanoseq$num_indels_per_cell_uci,col="darkgray")
		points(nanoseq[which(nanoseq$cell_type=="smuscleU"),"age"],nanoseq[which(nanoseq$cell_type=="smuscleU"),"num_indels_per_cell"],pch=20,col="darkolivegreen")
		points(nanoseq[which(nanoseq$cell_type=="smuscleC"),"age"],nanoseq[which(nanoseq$cell_type=="smuscleC"),"num_indels_per_cell"],pch=20,col="darkolivegreen2")
		abline(coef[1], coef[2], col="darkolivegreen",lty=2)		
		legend("topleft",legend=c("Bladder","Colon"),col=c("darkolivegreen","darkolivegreen2"),pch=15,cex=.6,bty = "n")
	} else {	
		all_ages <- data.frame(age = c(0:100))
		ci=predict(lm(num_indels_per_cell~age,data=nanoseq),all_ages,interval="confidence", level=0.95)
		polygon(x=c(all_ages$age, rev(all_ages$age)), y=c(ci[,3], rev(ci[,2])), col="grey90", border=NA)
		segments(nanoseq$age,nanoseq$num_indels_per_cell_lci,nanoseq$age,nanoseq$num_indels_per_cell_uci,col="darkgray")
		points(nanoseq[which(nanoseq$disease_state=="Healthy"),"age"],nanoseq[which(nanoseq$disease_state=="Healthy"),"num_indels_per_cell"],pch=20,col="firebrick")
		points(nanoseq[which(nanoseq$disease_state=="AD"),"age"],nanoseq[which(nanoseq$disease_state=="AD"),"num_indels_per_cell"],pch=20,col="firebrick2")
		coef = lm(nanoseq$num_indels_per_cell~nanoseq$age)$coefficients
		abline(coef[1], coef[2], col="firebrick",lty=2)
		legend("topleft",legend=c("Healthy","Alzheimer's disease"),col=c("firebrick","firebrick2"),pch=15,cex=.6,bty = "n")
	}
}


##########################################################################################
# Panels c and d (substitution & indel profiles) obtained with botseq_results_plotter.R 
# Code available at https://github.com/cancerit/NanoSeq
# and with SigProfilerMatrixGenerator

profiles = read.table("../DATA/profiles.tsv", sep="\t", header=T, row.names=1)
profiles = profiles[grep("neurons",rownames(profiles)),]
colnames(profiles) = gsub("\\.",">",colnames(profiles))
obs      = apply(profiles,2,sum)

colours        = rep(c("deepskyblue","black","firebrick2","gray","darkolivegreen3","rosybrown2"),each=16)
sub_vec        = c("C>A","C>G","C>T","T>A","T>C","T>G")
ctx_vec        = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
full_vec       = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
xstr           = paste(substr(full_vec,5,5), substr(full_vec,1,1), substr(full_vec,7,7), sep="")
ordered_names  = paste(xstr,">",rep(c("A","G","T","A","C","G"),each=16),sep="")

obs = obs[ordered_names]

par(mfrow=c(2,1))
par(mar=c(3,6,3,3))
par(mgp=c(3,1,0))
y = obs
maxy = max(y)
h = barplot(y, las=2, col=colours, border=NA, ylim=c(0,maxy*1.5), space=1, cex.names=0.6, names.arg=xstr, 
            ylab="Corrected mutation counts",main="Neurons",cex.main=.7)
for (j in 1:length(sub_vec)) {
 	xpos = h[c((j-1)*16+1,j*16)]
    rect(xpos[1]-0.5, maxy*1.2, xpos[2]+0.5, maxy*1.3, border=NA, col=colours[j*16])
 	text(x=mean(xpos), y=maxy*1.3, pos=3, label=sub_vec[j])
}    


##########################################################################################
# Panel e - Signature extraction

# NANOSEQ SIGNATURE ANALYSIS
# 2020/05/21

library(sigfit)
library(stringr)

# A) SIGNATURE EXTRACTION
###########################

# Load counts and discard urothelium, bladder and colibactin-exposed samples
profiles = read.table("../DATA/profiles.tsv", sep="\t", header=T, row.names=1)
counts = round(profiles[!grepl("(PD37449)|(urothel)|(bladder)|(botseq)", rownames(profiles)), ])
rownames(counts) = gsub("nanoseqv2", "nanoseq", gsub(" nanoseq", "_nanoseq", rownames(counts)))


# Merge counts by tissue and sequencing type
types = c("cordblood_nanoseq", "cordblood_stdseq", "grans_nanoseq",
          "(adultblood_stdseq)|(colonies_stdseq)|(grans_colonies)",
          "(coloniccryptS_nanoseq)|(coloniccryptM_nanoseq)", "coloniccrypts_stdseq",
          "smoothmuscleU_nanoseq", "smoothmuscleC_nanoseq",
          "(neuronsFC_nanoseq)|(neuronsFC-AD_nanoseq)")

counts.per.type = t(sapply(types, function(type) {
    colSums(counts[grepl(type, rownames(counts)), , drop=F])
}))

# Downsample to a maximum of 2000 mutations per catalogue
MAX.MUT = 2000
counts.per.type.down = t(apply(counts.per.type, 1, function(cnt) {
    if (sum(cnt) < MAX.MUT)
        cnt
    else
        round(cnt / sum(cnt) * MAX.MUT)
}))


# Extract 3 signatures from merged catalogues
fit.merged.3 = extract_signatures(counts.per.type.down, 3,
                                  iter=15000, warmup=5000, seed=0xC0FFEE)
# plot_all(fit.merged.3, out_path="SignatureFitting_Merged",
#          exp_cex_names=0.9, exp_margin_bottom=14)
signatures = retrieve_pars(fit.merged.3, "signatures")$mean
write.table(signatures,file="../DATA/3signatures.tsv",sep="\t",row.names=T,col.names=T,quote=F)

# B) SIGNATURE EXPOSURE ANALYSES
##################################

# Load sample info table, reload counts
nanoseq = read.table("../DATA/rates.tsv",header=T,sep="\t",row.names=1)
nanoseq = nanoseq[which(nanoseq$cell_type=="neurons"),]
profiles = read.table("../DATA/profiles.tsv", sep="\t", header=T, row.names=1)
profiles = profiles[grep("neurons",rownames(profiles)),]

# Fit the extracted signatures to each sample
fit.sample.3 = fit_signatures(profiles, signatures, iter=5000, seed=0xC0FFEE, cores=4)

plot_all(fit.sample.3, out_path="../DATA/SignatureFitting_AllNeuronSamples",
         exp_cex_names=0.9, exp_margin_bottom=14)


# Convert exposures to mutations per cell
# Burden per cell = num_muts / num_sites_called * genome_size
neurons.counts = rowSums(profiles[grep("neurons", rownames(profiles)), ])

neurons.burden = neurons.counts / nanoseq[names(neurons.counts),"total"] * genome_size
neurons.expos = retrieve_pars(fit.sample.3, "exposures")$mean[names(neurons.counts), ]
stopifnot(identical(names(neurons.burden), rownames(neurons.expos)))

neurons.burden.sigs = neurons.burden * neurons.expos
stopifnot(max(abs(rowSums(neurons.burden.sigs) - neurons.burden)) < 1e-10)

# Output results
save(neurons.expos, neurons.burden, neurons.burden.sigs,
     file="../DATA/Neurons_Exposures_Burden.RData")

cat("\nDone\n")

##########################################################################################
# Panel 3i - barplot exposure signatures in neurons
# Sort by age
# Then sort by healthy, AD
nanoseq = nanoseq[order(nanoseq$age),]
nanoseq = nanoseq[order(nanoseq$disease_state,decreasing=T),]
profiles = profiles[rownames(nanoseq),]
neurons.burden.sigs = neurons.burden.sigs[rownames(nanoseq),]
neurons.expos = neurons.expos[rownames(nanoseq),]
rownames(profiles) = paste(rownames(nanoseq),"-",nanoseq$age,"yo",sep="")
rownames(neurons.burden.sigs) = paste(rownames(nanoseq),"-",nanoseq$age,"yo",sep="")
rownames(neurons.expos) = paste(rownames(nanoseq),"-",nanoseq$age,"yo",sep="")

par(mar=c(10,4,4,4))
barplot(t(as.matrix(neurons.expos)),las=2,space=c(rep(.2,8),1,rep(.2,8)),
        ylab="Signature contribution",ylim=c(0,1.2),
        col=c("darkorange3","goldenrod4","darkseagreen4"))
legend("top",ncol=3,legend=c("SigA","SigB","SigC"),col=c("darkorange3","goldenrod4","darkseagreen4"),pch=15,cex=1,bty = "n")

##########################################################################################
# Panel 3n - barplot exposure signatures in smooth muscle
# Sort by age
# Then sort by colon/bladder
nanoseq = read.table("../DATA/rates.tsv",header=T,sep="\t",row.names=1)
nanoseq = nanoseq[grep("muscle",nanoseq$cell_type),]
profiles = read.table("../DATA/profiles.tsv", sep="\t", header=T, row.names=1)
profiles = profiles[grep("muscle",rownames(profiles)),]

nanoseq = nanoseq[order(nanoseq$age),]
nanoseq = nanoseq[order(nanoseq$cell_type,decreasing=T),]
profiles = profiles[rownames(nanoseq),]

fit.sample.3 = fit_signatures(profiles, signatures, iter=5000, seed=0xC0FFEE, cores=4)

# Convert exposures to mutations per cell
# Burden per cell = num_muts / num_sites_called * genome_size
smuscle.counts = rowSums(profiles[grep("muscle", rownames(profiles)), ])

smuscle.burden = smuscle.counts / nanoseq[names(smuscle.counts),"total"] * genome_size
smuscle.expos = retrieve_pars(fit.sample.3, "exposures")$mean[names(smuscle.counts), ]
stopifnot(identical(names(smuscle.burden), rownames(smuscle.expos)))

smuscle.burden.sigs = smuscle.burden * smuscle.expos
stopifnot(max(abs(rowSums(smuscle.burden.sigs) - smuscle.burden)) < 1e-10)

smuscle.burden.sigs = smuscle.burden.sigs[rownames(nanoseq),]
smuscle.expos = smuscle.expos[rownames(nanoseq),]
rownames(profiles) = paste(rownames(nanoseq),"-",nanoseq$age,"yo",sep="")
rownames(smuscle.burden.sigs) = paste(rownames(nanoseq),"-",nanoseq$age,"yo",sep="")
rownames(smuscle.expos) = paste(rownames(nanoseq),"-",nanoseq$age,"yo",sep="")

par(mar=c(10,4,4,4))
barplot(t(as.matrix(smuscle.expos)),las=2,space=c(rep(.2,7),1,rep(.2,4)),
        ylab="Signature contribution",ylim=c(0,1.2),
        col=c("darkorange3","goldenrod4","darkseagreen4"))
legend("top",ncol=3,legend=c("SigA","SigB","SigC"),col=c("darkorange3","goldenrod4","darkseagreen4"),pch=15,cex=1,bty = "n")


##########################################################################################
# Panel 3h, indel rates and gene expression
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Given dependence on heavy files and third-party data, the code below is not ready
# to run directly
# Files can be requested to fa8@sanger.ac.uk
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Trinucleotide frequencies and trinucleotide substitution counts for each gene were
# calculated as explained in the BURDEN_IN_SPECIFIC_REGIONS section of the code repository

# Gene expression values for neurons and smooth muscle taken from GTex:
# GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct
# (Muscle...Skeletal, Brain...Frontal.Cortex..BA9., respectively)

# Indel types annotated with SigProfileMatrixGenerator
# List of indel calls available in Suppl. Table 5
# DATA/Indels.neurons.annotated.tsv
# DATA/Indels.smooth_muscle.annotated.tsv
# Alternatively, indel types can be obtained as in Panel p

# The gene_prot_coords.tsv file was obtained from Ensembl/Biomart, and its format should
# be:
# gene    prot    chr     start   end     tstart  tend    strand
# ENSG00000223116         13      23551994        23552136        23551994        23552136        -1
# ENSG00000233440         13      23708313        23708703        23708313        23708703        1
# ENSG00000207157         13      23726725        23726825        23726725        23726825        -1
# ENSG00000229483         13      23743974        23744736        23743974        23744736        -1
# etc

library(dndscv)
library(Rsamtools)
library(GenomicRanges)
library(Biostrings)
library(epitools)

# Choose one below:
tipo = "SMOOTH_MUSCLE"
tipo = "NEURONS"

ONLY_LONG_INDELS = 0

if(tipo == "SMOOTH_MUSCLE") {
	setwd("../DATA/GENE_RATES/SMOOTH_MUSCLE")
} else if(tipo == "NEURONS") {
	setwd("../DATA/GENE_RATES/NEURONS")
}

files = list.files(".","gene_rates.trinuc.summary.tsv")

# Remove contaminated samples - for indels I didn't carry an in silico decontamination
files = setdiff(files, c("NP20.gene_rates.trinuc.summary.tsv","NP21.gene_rates.trinuc.summary.tsv","61.gene_rates.trinuc.summary.tsv"))

res = data.frame()
for(file in files) {
	cat(file,"...\n")
	j = read.table(file,sep="\t",header=T)
	if(nrow(res)==0) {
		res = j
	} else {
		res = rbind(res,j)
	}
}
colnames(res) = gsub("\\.",">",colnames(res))

all = res
all.bck = all


# Load expression
exp = read.table("../DATA/GENE_EXPRESSION/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct",sep="\t",header=T,quote="")
exp = exp[which(exp$Description %in% unique(all$gene)),]
exp2_smooth   = exp[,c("Name","Description","Muscle...Skeletal")]
exp2_neurons  = exp[,c("Name","Description","Brain...Frontal.Cortex..BA9.")]
colnames(exp2_smooth) = colnames(exp2_neurons) = c("Name","Description","expression")
if(tipo == "SMOOTH_MUSCLE") {
	exp = exp2_smooth
} else if(tipo == "NEURONS") {
	exp = exp2_neurons
}


# Set some things:
	nbins = 4
	exp_order = exp[order(exp$expression,decreasing=F),"Description"]
	
	tris = colnames(all)[8:39]
	subs = colnames(all)[40:135]
	cell_types = c("grans","islets","smooth","neurons")
	
	library(deconstructSigs)
	tri_genome = tri.counts.genome[,1]
	names(tri_genome) = rownames(tri.counts.genome)
	tri_genome = tri_genome[tris]
	tri_genome_freqs = tri_genome/sum(tri_genome)
	tri_prefix = sapply(subs,function(x) substr(x,1,3))


if(tipo == "NEURONS") {
	indels_info = read.table("../DATA/indels.neurons.annotated.tsv",     sep="\t",header=F)
} else if(tipo == "SMOOTH_MUSCLE") {
	indels_info = read.table("../DATA/indels.smooth_muscle.annotated.tsv",sep="\t",header=F)
}
colnames(indels_info) = c("file","chr","pos","info","ref","alt")
indels_info$len       = sapply(indels_info$info,function(x) unlist(strsplit(x,":"))[2])
indels_info$type      = sapply(indels_info$info,function(x) unlist(strsplit(x,":"))[3])
indels_info$fourth    = sapply(indels_info$info,function(x) unlist(strsplit(x,":"))[4])
indels_info$fifth     = sapply(indels_info$info,function(x) unlist(strsplit(x,":"))[5])

# To restrict the analysis to >1 bp indels
if(ONLY_LONG_INDELS==1) {
	indels_info = indels_info[which(indels_info$len>1),]
}

indels_info$indel_type = ""
indels_info[which(indels_info$type == "Del" & indels_info$len == 1),"indel_type"] = "Del1"
indels_info[which(indels_info$type == "Del" & indels_info$len >  1),"indel_type"] = "Del1+"
indels_info[which(indels_info$type == "Ins"                       ),"indel_type"] = "Ins"

##########################################################################################
# Load all the gene coordinates
	# Gene related data: lengths, gc content... and most importantly, which are protein-coding
	data("refcds_hg19", package="dndscv")
	gene_list = sapply(RefCDS, function(x) x$gene_name)
	protein_list = sapply(RefCDS, function(x) x$protein_id)
	gene2prot <- cbind(gene_list,protein_list)
	rownames(gene2prot) <- gene2prot[,1]
	gene2prot <- gene2prot[,2]
	cdss <- sapply(RefCDS,function(x) strsplit(as.character(x$seq_cds),split = "")[[1]])
	tmp_ <- sapply(cdss,table)
	gc   <- sapply(tmp_,function(x) (x["G"]+x["C"])/(x["G"]+x["C"]+x["A"]+x["T"]))
	names(gc) <- gene_list
	lens   <- sapply(cdss,length)
	
	# Also needed are: transcription strand, start & end of gene
	gene_names  <- sapply(RefCDS,function(x) x$gene_name )
	strands     <- sapply(RefCDS,function(x) x$strand    )
	gene_ids    <- sapply(RefCDS,function(x) x$gene_id   )
	prot_ids    <- sapply(RefCDS,function(x) x$protein_id)
	id2name     <- vector();	id2name[gene_ids] <- gene_names
	biomart     <- read.table("../DATA/gene_prot_coords.tsv",header=T,sep="\t")
	#biomart     <- read.table("/lustre/scratch119/casm/team268im/fa8/BOTSEQ/SUBST_BIASES_TRANSCRIPTION_ET_AL/gene_prot_coords.tsv",header=T,sep="\t")
	biomart     <- biomart[grep("ENSP",biomart$prot),]
	rownames(biomart) <- biomart$prot
	biomart$true_start <- biomart$start; biomart$true_end <- biomart$end;
	biomart[which(biomart$strand==-1),"true_start"] <- biomart[which(biomart$strand==-1),"tend"  ]
	biomart[which(biomart$strand==-1),"true_end"  ] <- biomart[which(biomart$strand==-1),"tstart"]
	tstarts <- biomart[prot_ids,"true_start"]
	tends   <- biomart[prot_ids,"true_end"  ]


	##########################################################################################
	# Entire gene regions +/-10Kb:
	biomart_prot2indexes <- c(1:nrow(biomart)); names(biomart_prot2indexes) <- biomart$prot
	
	genes_in_plus_indexes   <- which(strands== 1)
	chrs                    <- biomart[biomart_prot2indexes[prot_ids[genes_in_plus_indexes]],"chr"]
	tstarts                 <- biomart[biomart_prot2indexes[prot_ids[genes_in_plus_indexes]],"tstart"]
	tends                   <- biomart[biomart_prot2indexes[prot_ids[genes_in_plus_indexes]],"tend"]
	strands_                <- rep("+",length(chrs))
	genes_plus_tmp          <- data.frame(gene=gene_names[genes_in_plus_indexes],gene_id=gene_ids[genes_in_plus_indexes],chr=chrs,tstart=tstarts,tend=tends,strand=strands_,            type=rep("gene",length(genes_in_plus_indexes)))
	igr_in_plus_tmp_upst    <- data.frame(gene=gene_names[genes_in_plus_indexes],gene_id=gene_ids[genes_in_plus_indexes],chr=chrs,tstart=tstarts-10000-1,tend=tstarts-1,strand=strands_,type=rep("upst",length(genes_in_plus_indexes)))
	igr_in_plus_tmp_down    <- data.frame(gene=gene_names[genes_in_plus_indexes],gene_id=gene_ids[genes_in_plus_indexes],chr=chrs,tstart=tends+1,tend=tends+10000+1,strand=strands_,    type=rep("down",length(genes_in_plus_indexes)))
	
	genes_in_minus_indexes  <- which(strands==-1)
	chrs                    <- biomart[biomart_prot2indexes[prot_ids[genes_in_minus_indexes]],"chr"]
	tstarts                 <- biomart[biomart_prot2indexes[prot_ids[genes_in_minus_indexes]],"tstart"]
	tends                   <- biomart[biomart_prot2indexes[prot_ids[genes_in_minus_indexes]],"tend"]
	strands_                <- rep("-",length(chrs))
	genes_minus_tmp         <- data.frame(gene=gene_names[genes_in_minus_indexes],gene_id=gene_ids[genes_in_minus_indexes],chr=chrs,tstart=tstarts,tend=tends,strand=strands_,             type=rep("gene",length(genes_in_minus_indexes)))
	igr_in_minus_tmp_upst   <- data.frame(gene=gene_names[genes_in_minus_indexes],gene_id=gene_ids[genes_in_minus_indexes],chr=chrs,tstart=tends+1,tend=tends+10000+1,strand=strands_,     type=rep("upst",length(genes_in_minus_indexes)))
	igr_in_minus_tmp_down   <- data.frame(gene=gene_names[genes_in_minus_indexes],gene_id=gene_ids[genes_in_minus_indexes],chr=chrs,tstart=tstarts-10000-1,tend=tstarts-1,strand=strands_, type=rep("down",length(genes_in_minus_indexes)))
	
	whole_genes_df          <- rbind(genes_plus_tmp,genes_minus_tmp)
	whole_igr_df            <- rbind(igr_in_plus_tmp_upst,igr_in_plus_tmp_down,igr_in_minus_tmp_upst,igr_in_minus_tmp_down)
	
	whole_genes_df$id       <- paste(whole_genes_df$gene,whole_genes_df$type,sep=":")
	whole_igr_df$id         <- paste(whole_igr_df$gene,whole_igr_df$type,sep=":")
	rownames(whole_genes_df) <- whole_genes_df$id
	rownames(whole_igr_df  ) <- whole_igr_df$id


	# Rates by expression quartiles, showing types of indels 
	#  (Del1, Del1+, Ins)
	##########################################################################################
	cell_type = tipo
	all = all.bck[which(all.bck$type == "gene"),]
	##########################################################################################
	# Coverage by gene
	jorl = all[,8:(8+32-1)]
	genes = all$gene
	total_coverage = tapply(apply(jorl,1,sum),genes,sum)
	
	##########################################################################################
	# Intersect
	all_res = list()
	for(indel_type in c("Ins","Del1","Del1+")) {
		indels_gr         = GRanges(indels_info[which(indels_info$indel_type==indel_type),"chr"], IRanges(start=indels_info[which(indels_info$indel_type==indel_type),"pos"], end=indels_info[which(indels_info$indel_type==indel_type),"pos"]        ))
		genes_gr          = GRanges(whole_genes_df[,"chr"], IRanges(start=whole_genes_df[,"tstart"]-1, end=whole_genes_df[,"tend"]))
		indels_in_regions = as.data.frame(findOverlaps(indels_gr,genes_gr))
	
		# Counts per gene:
		kk = sort(table(whole_genes_df[indels_in_regions[,"subjectHits"],"gene"]))
		indel_rates = data.frame(cov=total_coverage,indels=0)
		rownames(indel_rates) = names(total_coverage)
		indel_rates[names(kk),"indels"]=kk

		##########################################################################################
		# Calculate values for each bin
		res = as.data.frame(matrix(nrow=0,ncol=6))
		colnames(res) = c("expr_bin","indels","sites","rate","rate_lci","rate_uci")
		for(nbin in c(1:nbins)) {
			start = round(length(exp_order)/nbins) * (nbin-1) + 1
			end   = start + round(length(exp_order)/nbins) 
			genes_in_bin = exp_order[start:end]
			
			res_ = vector()
			
			res_["expr_bin"]  = nbin
			res_["indels"]    = sum(indel_rates[which(rownames(indel_rates) %in% genes_in_bin),"indels"])
			res_["sites"]     = sum(indel_rates[which(rownames(indel_rates) %in% genes_in_bin),"cov"])
			res_["rate"]      = res_["indels"] / res_["sites"]
			res_["rate_lci"]  = poisson.test(res_["indels"])$conf.int[1] / res_["sites"]
			res_["rate_uci"]  = poisson.test(res_["indels"])$conf.int[2] / res_["sites"]
	
			if(nrow(res)==0) {
				res[1,] = res_
			} else {
				res = rbind(res,res_)
			}
		}
		all_res[[indel_type]] = res
	}
	
	rates = matrix(nrow=nbins,ncol=3)
	colnames(rates) = c("Ins","Del1","Del1+")
	ucis = lcis = rates
	for(bin in 1:nbins) {
		rates[bin,] = c(all_res[["Ins"]][bin,"rate"],all_res[["Del1"]][bin,"rate"],all_res[["Del1+"]][bin,"rate"])
		ucis[bin,] = c(all_res[["Ins"]][bin,"rate_uci"],all_res[["Del1"]][bin,"rate_uci"],all_res[["Del1+"]][bin,"rate_uci"])
		lcis[bin,] = c(all_res[["Ins"]][bin,"rate_lci"],all_res[["Del1"]][bin,"rate_lci"],all_res[["Del1+"]][bin,"rate_lci"])
	}
	plot(1:4, rep(max(ucis),4), type="n",xaxt="n",xlab="",ylim=c(0,max(ucis)),xlim=c(0.7,4.6),
	     ylab="Indel rate (cohort)",las=2,
	     frame.plot=F, main="Neurons")
	colours = vector(); colours[colnames(rates)] = c("darkolivegreen3","aquamarine3","aquamarine4")
	for(i in c(1:3)) {
		tipo = colnames(rates)[i]
		color = colours[tipo]
		xcoords = c(1:4)+(i*.1)
		segments(xcoords,lcis[,tipo],xcoords,ucis[,tipo],col="darkgray")
		points(xcoords,rates[,tipo],col=color,pch=20)
		segments(xcoords[1:3],rates[,tipo][1:3],xcoords[2:4],rates[,tipo][2:4],col=color)
	}
	legend("topleft",legend=colnames(rates),col=colours,pch=20,lwd=1,cex=.8)
	axis(side=1,at=c(1:4)+.2,labels = c("Q1","Q2","Q3","Q4"))


##########################################################################################
# Panels l and m: substitution and indel profiles
# Indel profile obtained with SigProfilerMatrixGenerator

profiles = read.table("../DATA/profiles.tsv", sep="\t", header=T, row.names=1)
profiles = profiles[grep("muscle",rownames(profiles)),]
colnames(profiles) = gsub("\\.",">",colnames(profiles))
obs      = apply(profiles,2,sum)

colours        = rep(c("deepskyblue","black","firebrick2","gray","darkolivegreen3","rosybrown2"),each=16)
sub_vec        = c("C>A","C>G","C>T","T>A","T>C","T>G")
ctx_vec        = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
full_vec       = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
xstr           = paste(substr(full_vec,5,5), substr(full_vec,1,1), substr(full_vec,7,7), sep="")
ordered_names  = paste(xstr,">",rep(c("A","G","T","A","C","G"),each=16),sep="")

obs = obs[ordered_names]

par(mfrow=c(2,1))
par(mar=c(3,6,3,3))
par(mgp=c(3,1,0))
y = obs
maxy = max(y)
h = barplot(y, las=2, col=colours, border=NA, ylim=c(0,maxy*1.5), space=1, cex.names=0.6, names.arg=xstr, 
            ylab="Corrected mutation counts",main="Smooth muscle",cex.main=.7)
for (j in 1:length(sub_vec)) {
 	xpos = h[c((j-1)*16+1,j*16)]
    rect(xpos[1]-0.5, maxy*1.2, xpos[2]+0.5, maxy*1.3, border=NA, col=colours[j*16])
 	text(x=mean(xpos), y=maxy*1.3, pos=3, label=sub_vec[j])
}    

##########################################################################################
# Panels o and p: substitution and indel rates across cell types

nanoseq = read.table("../DATA/rates.tsv",header=T,sep="\t",row.names=1)
nanoseq = nanoseq[setdiff(rownames(nanoseq),c("PD43976-Botseq","PD43976-Nanoseq")),]

j = nanoseq

j$mut_rate = j$burden
j$num_sites = j$total

j = j[which(j$age > 0),] # remove cord blood

j$donor = substr(rownames(j),1,7)
j$donor = paste(j$donor,j$cell_type,sep="-")
j = j[grep("PD37449",rownames(j),invert=T),                     ] # remove episodic colibactin
j = j[grep("PD40794_urothelium_nanoseqv2",rownames(j),invert=T),] # remove donor with SigA

# To collapse smuscle:
# j$donor = gsub("smuscle.","smuscle",j$donor);j$cell_type = gsub("smuscle.","smuscle",j$cell_type)

library(deconstructSigs)
genome_size = sum(tri.counts.genome)*2

j$muts_per_cell = j$mut_rate * genome_size
j$muts_per_year = j$muts_per_cell/j$age

# Do the linear regression to get confidence intervals
ordered_names = unique(j$cell_type)
result = matrix(nrow=length(ordered_names),ncol=3)
colnames(result) = c("muts_per_year","cilow","cihigh")
rownames(result) = ordered_names
for(cu in ordered_names) {
	rr = j[which(j$cell_type==cu),]
	cat("cu=",cu,"\n")
	print(rr[,c("muts_per_cell","age","donor")])
	if(length(unique(rr$donor)) == nrow(rr)) {
		model = lm(muts_per_cell ~ age-1, data=rr)
		ci = confint(model)
		result[cu,] = c(coef(model),as.vector(ci))
	} else {
		library(lme4)
		model = lmer(muts_per_cell ~ age + (age-1|donor) - 1,data=rr)
		ci = confint(model)["age",]
		result[cu,] = c(fixef(model),as.vector(ci))
	}
}
result[result<0] = 0
result = result[order(result[,1]),]

keep_order = rownames(result)

par(mfrow=c(1,2))
par(mgp=c(2, .7, 0))
bar = barplot(result[,1],names.arg=rownames(result),las=2,
              ylab="Number of substitutions per year", 
              ylim=c(0,max(result[,3])+max(result[,3])*.1),
              border=F,cex.names=.7,col="firebrick3")
segments(bar,result[,2],bar,result[,3],col="gray")

#########
# INDELS:
j = j[which(j$exclude_indels == "No"),]

########################
# Add indel information:
files = list.files("../DATA/INDEL_CALLS/",pattern="bulkmasked.tsv")
mat = data.frame()
for(file in files) {
	jjj = read.table(paste("../DATA/INDEL_CALLS/",file,sep=""),sep="\t",header=T)
	jjj = jjj[grep("PASS",jjj$FILTER),]
	sampleID = unlist(strsplit(file,"\\."))[1]
	sampleID = gsub("ARJ","PD43976",sampleID)
	jjj$sampleID = sampleID
	if(length(grep("(ARJ)|(grans)",sampleID)) != 0) {
		jjj$type = "grans"
	} else if(length(grep("neuron",sampleID)) != 0) {
		jjj$type = "neurons"
	} else if(length(grep("(colon)",sampleID)) != 0) {
		jjj$type = "colonic_crypt"
	} else if(length(grep("(muscle)",sampleID)) != 0) {
		jjj$type = "smooth_muscle"
	} else if(length(grep("(sperm)",sampleID)) != 0) {
		jjj$type = "sperm"
	} else if(length(grep("(urothel)",sampleID)) != 0) {
		jjj$type = "urothelium"
	} else {
		next
	}
	if(nrow(mat) == 0) {
		mat = jjj
	} else {
		mat = rbind(mat,jjj)
	}
}

# For multiallelic indels, I will consider only the first one
mat$to = sapply(mat$ALT,function(x) unlist(strsplit(x,","))[1])
mat$class = "DEL"
mat[which(nchar(mat$to)>nchar(mat$REF)),"class"] = "INS" 
mat[which(nchar(mat$to)+1<nchar(mat$REF)),"class"] = "LONG_DEL" 

# Excluding insertions
kk = t(table(mat[which(mat$class!="INS"),c("class","type")]))
kk / apply(kk,1,sum)

# Insertions vs deletions
mat$class2 = gsub("LONG_","",mat$class)
kk = t(table(mat[,c("class2","type")]))
kk / apply(kk,1,sum)

k = table(mat[,c("sampleID","class")])
k = k[rownames(j),]
j$DEL      = k[rownames(j),"DEL"]
j$LONG_DEL = k[rownames(j),"LONG_DEL"]
j$INS      = k[rownames(j),"INS"]

# To collapse smuscle:
# j$donor = gsub("smuscle.","smuscle",j$donor);j$cell_type = gsub("smuscle.","smuscle",j$cell_type)

j$DEL_per_cell       = j$DEL *      (genome_size/j$num_sites)
j$LONG_DEL_per_cell  = j$LONG_DEL * (genome_size/j$num_sites)
j$INS_per_cell       = j$INS *      (genome_size/j$num_sites)
j$DEL_per_year       = j$DEL *      (genome_size/j$num_sites) / j$age
j$LONG_DEL_per_year  = j$LONG_DEL * (genome_size/j$num_sites) / j$age
j$INS_per_year       = j$INS *      (genome_size/j$num_sites) / j$age

ordered_names = keep_order
result = matrix(nrow=length(ordered_names),ncol=3)
colnames(result) = c("INS","DEL","LONG_DEL")
rownames(result) = ordered_names
lowci  = result
hici   = result
for(cu in ordered_names) {
	for(tipo in colnames(result)) {
		rr = j[which(j$cell_type==cu),]
		field = paste(tipo,"_per_cell",sep="")
		
		if(length(unique(rr$donor)) == nrow(rr)) {
			model = lm(rr[,field] ~ rr[,"age"]-1, data=rr)
			ci = confint(model)
			result[cu,tipo] = coef(model)
			lowci[cu,tipo] = as.vector(ci)[1]
			hici[cu,tipo] = as.vector(ci)[2]
		} else {
			library(lme4)
			model = lmer(rr[,field] ~ age + (age-1|donor) - 1,data=rr)
			ci = confint(model)["age",]
			result[cu,tipo] = fixef(model)
			lowci[cu,tipo] = as.vector(ci)[1]
			hici[cu,tipo] = as.vector(ci)[2]
		}
	}
}
result[result<0] = 0
lowci[lowci<0] = 0

bar = barplot(t(result),
              beside=T,las=2,ylab="Indels per cell per year",
              col=c("darkolivegreen3","aquamarine3","aquamarine4"),
              border=NA,
              ylim=c(0,max(hici)))
legend("topleft",legend=c("Insertions","1 bp deletions", ">1 bp deletions"),col=c("darkolivegreen3","aquamarine3","aquamarine4"),pch=15)
segments(as.vector(bar),as.vector(t(lowci)),as.vector(bar),as.vector(t(hici)),col="gray")





