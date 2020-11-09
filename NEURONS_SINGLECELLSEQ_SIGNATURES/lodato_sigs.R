##########################################################################################
# This script processes Lodato et al's single cell sequencing, SNP-phased, mutation calls 
# (PMID:29217584) and extracts signatures with sigfit.
# This script also compares signatures to NanoSeq and to a signature specific of single 
# cell sequencing (Petljak et al (2019, PMID:30849372))
#
# Script by Yichen Wang and Federico Abascal
#
# Modify the genomeFile variable accordingly
##########################################################################################

library(readxl)
library(sigfit)
library("GenomicRanges")
library("Rsamtools")
library("MASS")
library(lsa) # for cosine similarity	
genomeFile = "/Users/fa8/Desktop/PANCANCER/TRANSCRIPTION_VS_INDELS_ANALYSIS/hs37d5.fa"

setwd("NEURONS_SINGLECELLSEQ_SIGNATURES")
# Supplementary tables from Lodato et al:
sheetnames <- excel_sheets("aao4426_TableS4_rev.xlsx")
mylist <- lapply(excel_sheets("aao4426_TableS4_rev.xlsx"), read_excel, path = "aao4426_TableS4_rev.xlsx")
names(mylist) <- sheetnames  # name the dataframes
list2env(mylist ,.GlobalEnv) # Bring the dataframes to the global environment

##########################################################################################
# Load all data and annotate mutations with trinucleotide/pyrimidine context (96 sub. types)
	by_donor=data.frame()
	for (cell in sheetnames){
	    sample_id=unique(mylist[[cell]][10])
	    mutations = cbind(mylist[[cell]][,1:4],mylist[[cell]][10])
	 
	    # 1. Annotating the trinucleotide context    
	    cat(nrow(mutations), "mutations loaded\n")
	    colnames(mutations) = c("chr","pos","ref","mut","sampleID")
	    cat("Annotating trinuc context...\n")
	    mutations = mutations[(mutations$ref %in% c("A","C","G","T")) & (mutations$mut %in% c("A","C","G","T")) & mutations$chr %in% c(1:22,"X","Y"),]
	    mutations$trinuc_ref = as.vector(scanFa(genomeFile, GRanges(mutations$chr, IRanges(mutations$pos-1, mutations$pos+1))))
	    
	    # 2. Annotating the mutation from the pyrimidine base
	    cat("Annotating pyr...\n")
	    ntcomp = c(T="A",G="C",C="G",A="T")
	    mutations$sub = paste(mutations$ref,mutations$mut,sep=">")
	    mutations$trinuc_ref_py = mutations$trinuc_ref
	    for (j in 1:nrow(mutations)) {
	        if (mutations$ref[j] %in% c("A","G")) { # Purine base
	            mutations$sub[j] = paste(ntcomp[mutations$ref[j]],ntcomp[mutations$mut[j]],sep=">")
	            mutations$trinuc_ref_py[j] = paste(ntcomp[rev(strsplit(mutations$trinuc_ref[j],split="")[[1]])],collapse="")
	        }
	    }
	    
	    
	    # 3. Counting subs
	    cat("Counting subs...\n")
	    freqs = table(paste(mutations$sub,paste(substr(mutations$trinuc_ref_py,1,1),substr(mutations$trinuc_ref_py,3,3),sep="-"),sep=","))
	    sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
	    ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
	    full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
	    freqs_full = freqs[full_vec]; freqs_full[is.na(freqs_full)]=0 ; names(freqs_full) = full_vec
	    
	    xstr = paste(substr(full_vec,5,5), substr(full_vec,1,1), substr(full_vec,7,7), sep="")
		xstr2 = paste(xstr,">",substr(full_vec,3,3),sep="")
	    
	    cell_by_donor<-t(as.data.frame(freqs_full))
	    cell_by_donor[2,]=as.numeric(cell_by_donor[2,])
	    colnames(cell_by_donor)=cell_by_donor[1,]
	    by_donor=rbind(by_donor,cell_by_donor)
	}
	by_donor = by_donor[grep("Var",rownames(by_donor),invert=T),]
	colnames(by_donor) = xstr2
	for(j in 1:96) {	by_donor[,j] = as.numeric(by_donor[,j])	}


##########################################################################################
# Plot the entire profile
	pdf("Lodato.spectrum_cosmiclike.pdf",width=12,height=4)
	colvec = rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16)
	y = freqs_full; maxy = max(y)
	h = barplot(y, las=2, col=colvec, border=NA, ylim=c(0,maxy*1.5), space=1, cex.names=0.6, names.arg=xstr, ylab="Number mutations", main="Lodato et al, all neurons")
	for (j in 1:length(sub_vec)) {
	    xpos = h[c((j-1)*16+1,j*16)]
	    rect(xpos[1]-0.5, maxy*1.2, xpos[2]+0.5, maxy*1.3, border=NA, col=colvec[j*16])
	    text(x=mean(xpos), y=maxy*1.3, pos=3, label=sub_vec[j])
	}    
	dev.off()

##########################################################################################
# Extract signatures with sigfit
	library(sigfit)
	data("cosmic_signatures_v3")	
	profiles = data.matrix(by_donor)
	### Select number of components:
	# nopriors_multinomial_guess_numcpts = extract_signatures(profiles,nsignatures=2:10,seed=1469,iter=1000)	
	# plot_gof(nopriors_multinomial_guess_numcpts, stat = "cosine")
	# plot_gof(nopriors_multinomial_guess_numcpts, stat = "L2")
	three_sig_extraction = extract_signatures(profiles,nsignatures=3,seed=1469,iter=10000)
	signatures3 = retrieve_pars(three_sig_extraction,par="signatures")
	plot_spectrum(signatures3,pdf_path="3cpts_signatures.LODATO.pdf")
	### Exposures:
	# mcmc_samples_refit <- fit_signatures(counts = profiles,
	#                                      signatures = signatures3,
	#                                      iter = 2000,
	#                                      warmup = 1000)
	# exposures <- retrieve_pars(mcmc_samples_refit, "exposures")
	# plot_exposures(counts = profiles, exposures = exposures,
	#                pdf_path = "Exposures3SIGS-LODATO.pdf")
	
##########################################################################################
# Comparison with NanoSeq data:
	profiles = read.table("triprofiles.tsv",sep="\t",header=T)
	profiles = profiles[grep("neuron",rownames(profiles)),]
	cosine(apply(profiles,2,sum),as.vector(as.matrix(signatures3$mean[1,]))) # 0.6593344
	cosine(apply(profiles,2,sum),as.vector(as.matrix(signatures3$mean[2,]))) # 0.5351618
	cosine(apply(profiles,2,sum),as.vector(as.matrix(signatures3$mean[3,]))) # 0.9505185 **


##########################################################################################
# Comparison with single-cell specific signatures from Petljak et al (2019, PMID:30849372):
	mias = read.table("mia_signatures.tsv",sep="\t",header=T);
	mia_scf = mias[,"SBS.sc_F"]
	names(mia_scf) = paste(mias$Mutation.Subtype,substr(mias$Mutation.Type,3,3),sep=">")
	cosine(mia_scf,as.vector(as.matrix(signatures3$mean[1,]))) # 0.412436
	cosine(mia_scf,as.vector(as.matrix(signatures3$mean[2,]))) # 0.965008 **
	cosine(mia_scf,as.vector(as.matrix(signatures3$mean[3,]))) # 0.427658
	
	pdf("scF.spectrum_cosmiclike.pdf",width=12,height=4)
	colvec = rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16)
	y = mia_scf; maxy = max(y)
	h = barplot(y, las=2, col=colvec, border=NA, ylim=c(0,maxy*1.5), space=1, cex.names=0.6, names.arg=xstr, ylab="Number mutations",main="scF (single cell seq specific signature)")
	for (j in 1:length(sub_vec)) {
	    xpos = h[c((j-1)*16+1,j*16)]
	    rect(xpos[1]-0.5, maxy*1.2, xpos[2]+0.5, maxy*1.3, border=NA, col=colvec[j*16])
	    text(x=mean(xpos), y=maxy*1.3, pos=3, label=sub_vec[j])
	}    
	dev.off()
	
##########################################################################################
# Reload Lodato et al data by neuron instead of by donor:
	by_neuron=data.frame()
	for (cell in sheetnames){
	    sample_ids=as.vector(as.matrix(unique(mylist[[cell]][11])))
	    mutations = cbind(mylist[[cell]][,1:4],mylist[[cell]][11])
	   
	    # 1. Annotating the trinucleotide context
	    cat(nrow(mutations), "mutations loaded\n")
	    colnames(mutations) = c("chr","pos","ref","mut","sampleID")
	    cat("Annotating trinuc context...\n")
	    mutations = mutations[(mutations$ref %in% c("A","C","G","T")) & (mutations$mut %in% c("A","C","G","T")) & mutations$chr %in% c(1:22,"X","Y"),]
	    mutations$trinuc_ref = as.vector(scanFa(genomeFile, GRanges(mutations$chr, IRanges(mutations$pos-1, mutations$pos+1))))	    
	    
	    # 2. Annotating the mutation from the pyrimidine base
	    cat("Annotating pyr...\n")
	    ntcomp = c(T="A",G="C",C="G",A="T")
	    mutations$sub = paste(mutations$ref,mutations$mut,sep=">")
	    mutations$trinuc_ref_py = mutations$trinuc_ref
	    for (j in 1:nrow(mutations)) {
	        if (mutations$ref[j] %in% c("A","G")) { # Purine base
	            mutations$sub[j] = paste(ntcomp[mutations$ref[j]],ntcomp[mutations$mut[j]],sep=">")
	            mutations$trinuc_ref_py[j] = paste(ntcomp[rev(strsplit(mutations$trinuc_ref[j],split="")[[1]])],collapse="")
	        }
	    }
	    
	    # 3. Counting subs
	    cat("Counting subs...\n")
	    #freqs = table(paste(mutations$sub,paste(substr(mutations$trinuc_ref_py,1,1),substr(mutations$trinuc_ref_py,3,3),sep="-"),sep=","))
	    for(ss in sample_ids) {
	    	submuts = mutations[which(mutations$sampleID %in% ss),]
	    	freqs = table(paste(submuts$sub,paste(substr(submuts$trinuc_ref_py,1,1),substr(submuts$trinuc_ref_py,3,3),sep="-"),sep=","))
	    	sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
	    	ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
	    	full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
	    	freqs_full = freqs[full_vec]; freqs_full[is.na(freqs_full)]=0 ; names(freqs_full) = full_vec
	    	
	    	xstr = paste(substr(full_vec,5,5), substr(full_vec,1,1), substr(full_vec,7,7), sep="")
			xstr2 = paste(xstr,">",substr(full_vec,3,3),sep="")
	    	
	    	cell_by_neuron<-t(as.data.frame(freqs_full))
	    	cell_by_neuron[2,]=as.numeric(cell_by_neuron[2,])
	    	colnames(cell_by_neuron)=cell_by_neuron[1,]
	    	rownames(cell_by_neuron) = c("KK",ss)
	    	by_neuron=rbind(by_neuron,cell_by_neuron)
	    }
	}
	
	by_neuron = by_neuron[grep("KK",rownames(by_neuron),invert=T),]
	colnames(by_neuron) = xstr2
	for(j in 1:96) {	by_neuron[,j] = as.numeric(by_neuron[,j])	}



##########################################################################################
# Calculate exposures to LDA, LDB and LDC for each neuron
	profiles = data.matrix(by_neuron)
	mcmc_samples_refit <- fit_signatures(counts = profiles,
	                                     signatures = signatures3,
	                                     iter = 2000,
	                                     warmup = 1000)
	exposures <- retrieve_pars(mcmc_samples_refit, "exposures")

	exposure_means = exposures$mean
	exposure_means$donor = sapply(rownames(exposure_means),function(x) unlist(strsplit(x,"_"))[1])
	colnames(exposure_means) = c("SignatureC","SintaureB","SignatureA","donor")
	donors = unique(exposure_means$donor)
	library(RColorBrewer)
	colors = vector()
	colors[donors] = terrain.colors(length(donors))
	counts_per_sig = exposure_means[,1:3] * apply(profiles,1,sum)
	counts_per_sig$donor = exposure_means$donor
		
	# Plot exposures
	par(mfrow=c(1,1))
	box = boxplot(exposure_means[,c(3,2,1)],notch=T,outline=F,at=c(2,3,4),ylab="Signature exposure")
	points(rep(3,nrow(exposure_means))+jitter(rep(1,nrow(exposure_means)),10),exposure_means[,1],col=rgb(.1,.1,.1,.2),pch=20)
	points(rep(2,nrow(exposure_means))+jitter(rep(1,nrow(exposure_means)),10),exposure_means[,2],col=rgb(.1,.1,.1,.2),pch=20)
	points(rep(1,nrow(exposure_means))+jitter(rep(1,nrow(exposure_means)),10),exposure_means[,3],col=rgb(.1,.1,.1,.2),pch=20)


##########################################################################################








