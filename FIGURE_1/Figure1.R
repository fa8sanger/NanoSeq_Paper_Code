##########################################################################################
# Panel b - Burden comparisons - boxplot colonies + botseq + nanoseq

library(deconstructSigs)
genome_size = sum(tri.counts.genome)*2

setwd("~/Desktop/BOTSEQ/PAPER/NEW_CODE_AND_MATRICES/FIGURE_1")
nanoseq                       = read.table("../DATA/rates.tsv",sep="\t",header=T)
nanoseq$num_muts_per_cell     = nanoseq$burden * genome_size
nanoseq$num_muts_per_cell_lci = nanoseq$burden_lci * genome_size
nanoseq$num_muts_per_cell_uci = nanoseq$burden_uci * genome_size

rr = read.table("../DATA/rates_caveman.tsv",sep="\t",header=F)
colnames(rr) = c("kk","file","sample","muts_caveman","c20_genome","muts","corr_muts","c20_nanoseq_genome","obs_rate","corr_rate")
files = unique(rr$file)
ages = vector()
ages[files] = c(59,63,54,38,36,0,29,38,48,81,0,77)
res = rr[which(rr$file == file),]
res = res[which(res$c20_nanoseq_genome>200e6),]
cat(median(res$c20_genome),"\n")
res$muts_per_cell = res$corr_rate * genome_size
cat(file,":",ages[file],"\n")

file = files[1]    
par(mar=c(6,4,3,3))
par(mgp=c(1.5,.7,0))
plot(NULL,NULL,xlim=c(0,3.5),ylim=c(0,nanoseq[which(nanoseq$sampleID=="PD43976-Botseq"),"num_muts_per_cell_uci"]),
     ylab="Number of mutations per cell",cex=.7,
     main="Grans ARJ",xaxt='n',xlab="",bty="n",frame.plot=F,cex.main=.6,cex.lab=.6,cex.axis=.6)

	boxplot(res$muts_per_cell,notch=T,
    	    cex.lab=.7,cex.names=.7,cex.main=.7,cex.axis=.7,
    	    add=T,at=0.75,frame.plot=F,axes=F,width=1.5,bty="n",col="darkcyan",
    	    pch=20,outcol="gray")
	rect(1.75,0,2.25,nanoseq[which(nanoseq$sampleID=="PD43976-Botseq"),"num_muts_per_cell"],
		col="firebrick", border=F)
	segments(2.0,nanoseq[which(nanoseq$sampleID=="PD43976-Botseq"),"num_muts_per_cell_lci"],
	         2.0,nanoseq[which(nanoseq$sampleID=="PD43976-Botseq"),"num_muts_per_cell_uci"],col="black",lwd=2)
	rect(2.75,0,3.25,nanoseq[which(nanoseq$sampleID=="PD43976-Nanoseq"),"num_muts_per_cell"],
		col="aquamarine4", border=F)
	segments(3.0,nanoseq[which(nanoseq$sampleID=="PD43976-Nanoseq"),"num_muts_per_cell_lci"],
	         3.0,nanoseq[which(nanoseq$sampleID=="PD43976-Nanoseq"),"num_muts_per_cell_uci"],col="black",lwd=2)
	axis(side=1, at=c(0.75,2,3), labels=c("Colonies","BotSeq","NanoSeq"), 
	     pos=0,las=1,cex=.6,cex.axis=.6,xlim=c(0,3.5))
	segments(0,0,3.5,0,col="black")


##################################################################################################
# Panel c - Trinucleotide substitution profiles: colonies, botseq, nanoseq
library(lsa)

profiles     = read.table("../DATA/profiles.tsv",sep="\t",header=T,row.names=1)
col_names    = colnames(profiles)
nanoseq_tri  = as.vector(as.matrix(profiles["PD43976_grans_nanoseq", ]))
colonies_tri = as.vector(as.matrix(profiles["PD43976_grans_colonies",]))
botseq_tri   = as.vector(as.matrix(profiles["PD43976_grans_botseq",  ]))
names(nanoseq_tri) = names(colonies_tri) = names(botseq_tri) = gsub("\\.",">",col_names)

	
	obs_botseq  = botseq_tri
	obs_nanoseq = nanoseq_tri
	ref         = colonies_tri

	# Calculate cosine similarity
	obs_cosine_sim_botseq  = cosine(obs_botseq, ref)
	obs_cosine_sim_nanoseq = cosine(obs_nanoseq,ref)
	
	# Calculate expected cosine similarity
	# Botseq
	cosines <- vector()
	simuls = list()
	SIMUL_SIZE = 1000
	for(i in c(1:SIMUL_SIZE)) {
		# Simulate dataset of size round(sum(obs))
		size = round(sum(obs_botseq))
		simul_ = sample(names(ref),size=size,p=ref/sum(ref),replace=T)
		simul  = table(simul_)
		simul[setdiff(names(ref),names(simul))] = 0
		simul  = as.vector(simul[names(ref)])
		names(simul) = names(obs_botseq)
		simuls[[i]] = simul
		#ref_tmp = ref - simul
		#cosines[i] = cosine(simul,ref_tmp)
		cosines[i] = cosine(simul,ref)
	}
	expected_cosine_sim_botseq       = mean(cosines) 
	stdev_expected_cosine_sim_botseq = sd(cosines)  
	lci_95_botseq                    = sort(cosines)[round(SIMUL_SIZE*0.025)] 
	uci_95_botseq                    = sort(cosines)[round(SIMUL_SIZE*0.975)] 
	cosines_botseq                   = cosines
	
	# Nanoseq
	cosines <- vector()
	simuls = list()
	SIMUL_SIZE = 1000
	for(i in c(1:SIMUL_SIZE)) {
		# Simulate dataset of size round(sum(obs))
		size = round(sum(obs_nanoseq))
		simul_ = sample(names(ref),size=size,p=ref/sum(ref),replace=T)
		simul  = table(simul_)
		simul[setdiff(names(ref),names(simul))] = 0
		simul  = as.vector(simul[names(ref)])
		names(simul) = names(obs_botseq)
		simuls[[i]] = simul
		#ref_tmp = ref - simul
		#cosines[i] = cosine(simul,ref_tmp)
		cosines[i] = cosine(simul,ref)
	}
	expected_cosine_sim_nanoseq       = mean(cosines) 
	stdev_expected_cosine_sim_nanoseq = sd(cosines)  
	lci_95_nanoseq                    = sort(cosines)[round(SIMUL_SIZE*0.025)] 
	uci_95_nanoseq                    = sort(cosines)[round(SIMUL_SIZE*0.975)] 
	cosines_nanoseq                   = cosines
	
	
	colours        = rep(c("deepskyblue","black","firebrick2","gray","darkolivegreen3","rosybrown2"),each=16)
	colours        = rep(c( rgb(34/255,159/255,198/255), rgb(26/255,26/255,26/255),
                            rgb(201/255,93/255,94/255),  rgb(178/255,182/255,180/255),
                            rgb(153/255,208/255,62/255), rgb(217/255,190/255,217/255)), each=16)

	sub_vec        = c("C>A","C>G","C>T","T>A","T>C","T>G")
	ctx_vec        = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
	full_vec       = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
	xstr           = paste(substr(full_vec,5,5), substr(full_vec,1,1), substr(full_vec,7,7), sep="")
	ordered_names  = paste(xstr,">",rep(c("A","G","T","A","C","G"),each=16),sep="")
	obs_botseq     = obs_botseq[ordered_names]
	obs_nanoseq    = obs_nanoseq[ordered_names]
	ref            = ref[ordered_names]
	
	######################################################################################
	# And plot it
	dev.new(width=8,height=7)
	par(mfrow=c(3,1))
	par(mar=c(2,4,5,3))
	par(mgp=c(3,1,0)) # or 3,1,0?
	
	#pdf("cord_blood.96_profile-profile_sigs.cosine.pdf",width=10,height=6)
	y = ref
	maxy = max(y)
	names(y) = NULL
	h = barplot(y, las=2, col=colours, border=NA, ylim=c(0,maxy*1.5), space=0.1, cex.names=0.6, names.arg=NULL, #names.arg=xstr, 
	            ylab="Mutation counts",main="Colonies")
	for (j in 1:length(sub_vec)) {
	 	xpos = h[c((j-1)*16+1,j*16)]
	    rect(xpos[1]-0.5, maxy*1.2, xpos[2]+0.5, maxy*1.3, border=NA, col=colours[j*16])
	 	text(x=mean(xpos), y=maxy*1.3, pos=3, label=sub_vec[j])
	}    
	y = obs_botseq
	names(y) = NULL
	maxy = max(y)
	h = barplot(y, las=2, col=colours, border=NA, ylim=c(0,maxy*1.5), space=.1, cex.names=0.6, names.arg=NULL, #names.arg=xstr, 
	            ylab="Corr. mutation counts",main=paste("Botseq","\nobs_cos_sim=",round(obs_cosine_sim_botseq,2),"\nexp_cos_sim=",round(expected_cosine_sim_botseq,2)," [CI95: ",round(lci_95_botseq,2),"-",round(uci_95_botseq,2),"]",sep=""),cex.main=.7)
	for (j in 1:length(sub_vec)) {
	 	xpos = h[c((j-1)*16+1,j*16)]
	    #rect(xpos[1]-0.5, maxy*1.2, xpos[2]+0.5, maxy*1.3, border=NA, col=colours[j*16])
	 	#text(x=mean(xpos), y=maxy*1.3, pos=3, label=sub_vec[j])
	}    
	y = obs_nanoseq
	maxy = max(y)
	h = barplot(y, las=2, col=colours, border=NA, ylim=c(0,maxy*1.5), space=.1, cex.names=0.6, names.arg=xstr, 
	            ylab="Corr. mutation counts",main=paste("Nanoseq","\nobs_cos_sim=",round(obs_cosine_sim_nanoseq,2),"\nexp_cos_sim=",round(expected_cosine_sim_nanoseq,2)," [CI95: ",round(lci_95_nanoseq,2),"-",round(uci_95_nanoseq,2),"]",sep=""),cex.main=.7)
	for (j in 1:length(sub_vec)) {
	 	xpos = h[c((j-1)*16+1,j*16)]
	    #rect(xpos[1]-0.5, maxy*1.2, xpos[2]+0.5, maxy*1.3, border=NA, col=colours[j*16])
	 	#text(x=mean(xpos), y=maxy*1.3, pos=3, label=sub_vec[j])
	}    
	
## boxplot(list(botseq=cosines_botseq,nanoseq=cosines_nanoseq),col="darkcyan",ylab="Cosine similarity",ylim=c(0,1))
## points(1,obs_cosine_sim_botseq,col="firebrick",pch=20,cex=2)
## points(2,obs_cosine_sim_nanoseq,col="aquamarine4",pch=20,cex=2)


##################################################################################################
# Fig 1d: asymmetries - plots obtained with botseq_results_plotter.R 
# Code available at https://github.com/cancerit/NanoSeq


##################################################################################################
# Fig 1f: sperm rates
# Substitutions per year obtained from the literature
	exp1 = 2.01
	lci1 = 1.68
	uci1 = 2.34
	exp2 = 2.87
	lci2 = 2.11
	uci2 = 3.64

library(deconstructSigs)
genome_size = sum(tri.counts.genome)*2

# Rate plot exluding the 73 yo
nanoseq               = read.table("../DATA/rates.tsv",sep="\t",header=T)
nanoseq               = nanoseq[which(nanoseq$cell_type=="sperm"&nanoseq$age==21),]
nanoseq$muts_per_cell = nanoseq$burden     * genome_size / 2 # sperm cells are haploid
nanoseq$lci_per_cell  = nanoseq$burden_lci * genome_size / 2 
nanoseq$uci_per_cell  = nanoseq$burden_uci * genome_size / 2

plot(NULL,NULL,xlim=c(0,13),ylim=c(0,max(max(nanoseq$uci_per_cell),3.64*21)),bty="n",frame.plot=F,xaxt="n",
     xlab="",ylab="Mutations per cell",main="Sperm")
rect(0.5,lci2*21,1.5,uci2*21,col=adjustcolor( "darkgoldenrod", alpha.f = 0.5),border=F)
rect(2.0,lci1*21,3.0,uci1*21,col=adjustcolor( "cyan4",         alpha.f = 0.5),border=F)
segments(0.5,exp2*21,1.5,exp2*21,col="darkgoldenrod",lwd=2)
segments(2.0,exp1*21,3.0,exp1*21,col="cyan4",lwd=2)

for(i in c(1:nrow(nanoseq))) {
	rect(4.6+i-1,
	     nanoseq[i,"muts_per_cell"],
	     5.4+i-1,
	     0,
	     col="aquamarine4",
	     border=F)
}
segments(c(5:(5+nrow(nanoseq)-1)),nanoseq$lci_per_cell,c(5:(5+nrow(nanoseq)-1)),nanoseq$uci_per_cell,lwd=1,col="black")

axis(side=1, at=c(1,2.5,c(5:(5+nrow(nanoseq)-1))), labels=c("Kong 2012","Rahbari 2016",rep("Nanoseq-21yo",nrow(nanoseq))), 
     pos=c(1,2.5,c(5:(5+nrow(nanoseq)-1))),las=2)





##################################################################################################
# Fig 1g: cord blood burdens
library(deconstructSigs)
genome_size = sum(tri.counts.genome)*2

# Colonies rates:
rr = read.table("../DATA/rates_caveman.tsv",sep="\t",header=F)
colnames(rr) = c("kk","file","sample","muts_caveman","c20_genome","muts","corr_muts","c20_nanoseq_genome","obs_rate","corr_rate")
files = unique(rr$file)
ages = vector()
ages[files] = c(59,63,54,38,36,0,29,38,48,81,0,77)
files = files[c(6,11)]
all_num_muts = vector()
corr_muts    = vector()
obs_muts     = vector()
c20_nanoseq_genomes = vector()
	for(file in files) {
	    res = rr[which(rr$file == file),]
	    res = res[which(res$c20_nanoseq_genome>200e6),]
	    cat(median(res$c20_genome),"\n")
		# boxplot(res[,c("muts","muts_per_cell")],main=file,
		res$muts_per_cell = res$corr_rate * genome_size
		cat(file,":",ages[file],"\n")
		all_num_muts = c(all_num_muts,res$muts_per_cell)
		corr_muts    = c(corr_muts,res$corr_muts)
		obs_muts     = c(obs_muts,res$muts)
		c20_nanoseq_genomes = c(c20_nanoseq_genomes,res$c20_nanoseq_genome)
	}
rate = sum(corr_muts)/sum(c20_nanoseq_genomes)
cis  = poisson.test(sum(obs_muts))$conf.int[c(1,2)]/sum(c20_nanoseq_genomes) * sum(corr_muts)/sum(obs_muts)
mean_muts = rate * genome_size
ci_muts   = cis * genome_size

# NanoSeq rates:
nanoseq                = read.table("../DATA/rates.tsv",sep="\t",header=T)
nanoseq                = nanoseq[which(nanoseq$age==0),]
nanoseq_correction     = 1/(1-0.24) # 24% percentage of mutations lost. I need to multiply by nanoseq_correction
nanoseq$mut_rate_final = nanoseq$burden        * nanoseq_correction
nanoseq$lci_final      = nanoseq$burden_lci    * nanoseq_correction
nanoseq$uci_final      = nanoseq$burden_uci    * nanoseq_correction
nanoseq$muts_per_cell  = nanoseq$mut_rate_final* genome_size
nanoseq$lci_final      = nanoseq$lci_final     * genome_size
nanoseq$uci_final      = nanoseq$uci_final     * genome_size

par(mar=c(10,4,3,3))
plot(NULL,NULL,xlim=c(0,10),ylim=c(0,max(nanoseq$uci_final)),
     ylab="Number of mutations per cell",cex=.7,
     main="Cord blood",xaxt='n',xlab="",bty="n",frame.plot=F)
boxplot(all_num_muts,at = 1,add = T, col="darkcyan",width=1.5,notch=T,bty="n",frame.plot=F,
         pch=20,outcol="gray")
segments(1,ci_muts[1],1,ci_muts[2],col="red")
points(1,mean_muts,pch=20,col="red",cex=.7)

colours = vector()
colours[1:6] = "aquamarine4"
colours[7]   = "aquamarine3"
for(i in c(1:nrow(nanoseq))) {
	rect(2.6+i-1,
	     nanoseq[i,"burden"] * genome_size,
	     3.4+i-1,
	     0,
	     col=colours[i],
	     border=F)
	rect(2.6+i-1,
	     nanoseq[i,"muts_per_cell"],
	     3.4+i-1,
	     nanoseq[i,"burden"] * genome_size,
	     col="azure3",
	     border=F)
}
segments(c(3:(3+nrow(nanoseq)-1)),nanoseq$lci_final,c(3:(3+nrow(nanoseq)-1)),nanoseq$uci_final,lwd=1,col="black")

axis(side=1, at=c(1,c(3:(3+nrow(nanoseq)-1))), labels=c("Colonies",c(rep("Nanoseq-S1",6),"Nanoseq-S2")), 
     pos=c(1,c(3:(3+nrow(nanoseq)-1))),las=2)


##################################################################################################
# Fig 1h: cord blood profiles

profiles     = read.table("../DATA/profiles.tsv",sep="\t",header=T,row.names=1)
col_names    = colnames(profiles)
colonies_prof1 = as.vector(as.matrix(profiles["CB001_cordblood_stdseq",]))
colonies_prof2 = as.vector(as.matrix(profiles["CB002_cordblood_stdseq",]))
nanoseq_prof1 = as.vector(as.matrix(profiles["XN3325_cordblood_nanoseq",]))
nanoseq_prof2 = as.vector(as.matrix(profiles["XN3326_cordblood_nanoseq",]))
names(colonies_prof1) = names(colonies_prof1) = names(nanoseq_prof1) = names(nanoseq_prof2) = gsub("\\.",">",col_names)

colonies_prof = colonies_prof1 + colonies_prof2
nanoseq_prof  = nanoseq_prof1  + nanoseq_prof2 

library(lsa)

obs  = nanoseq_prof
ref  = colonies_prof
obs = obs[names(ref)]

# Calculate cosine similarity
obs_cosine_sim  = cosine(obs, ref)

# Calculate expected cosine similarity
cosines <- vector()
simuls = list()
SIMUL_SIZE = 1000
for(i in c(1:SIMUL_SIZE)) {
	# Simulate dataset of size round(sum(obs))
	size = round(sum(obs))
	simul_ = sample(names(ref),size=size,p=ref/sum(ref),replace=T)
	simul  = table(simul_)
	simul[setdiff(names(ref),names(simul))] = 0
	simul  = as.vector(simul[names(ref)])
	names(simul) = names(obs)
	simuls[[i]] = simul
	#ref_tmp = ref - simul
	#cosines[i] = cosine(simul,ref_tmp)
	cosines[i] = cosine(simul,ref)
}
expected_cosine_sim       = mean(cosines) 
stdev_expected_cosine = sd(cosines)  
lci_95                    = sort(cosines)[round(SIMUL_SIZE*0.025)] 
uci_95                    = sort(cosines)[round(SIMUL_SIZE*0.975)] 
cosines                   = cosines

colours        = rep(c("deepskyblue","black","firebrick2","gray","darkolivegreen3","rosybrown2"),each=16)
colours        = rep(c( rgb(34/255,159/255,198/255), rgb(26/255,26/255,26/255),
                        rgb(201/255,93/255,94/255),  rgb(178/255,182/255,180/255),
                        rgb(153/255,208/255,62/255), rgb(217/255,190/255,217/255)), each=16)

sub_vec        = c("C>A","C>G","C>T","T>A","T>C","T>G")
ctx_vec        = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
full_vec       = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
xstr           = paste(substr(full_vec,5,5), substr(full_vec,1,1), substr(full_vec,7,7), sep="")
ordered_names  = paste(xstr,">",rep(c("A","G","T","A","C","G"),each=16),sep="")
obs            = obs[ordered_names]
ref            = ref[ordered_names]

# And plot it
dev.new(width=5,height=3)
par(mfrow=c(2,1))
par(mar=c(2,4,5,3))
par(mgp=c(3,1,0)) # or 3,1,0?

y = ref
maxy = max(y)
names(y) = NULL
h = barplot(y, las=2, col=colours, border=NA, ylim=c(0,maxy*1.5), space=0.1, cex.names=0.6, names.arg=NULL, #names.arg=xstr, 
            ylab="Mutation counts",main="Colonies")
for (j in 1:length(sub_vec)) {
 	xpos = h[c((j-1)*16+1,j*16)]
    rect(xpos[1]-0.5, maxy*1.2, xpos[2]+0.5, maxy*1.3, border=NA, col=colours[j*16])
 	text(x=mean(xpos), y=maxy*1.3, pos=3, label=sub_vec[j])
}    
y = obs
maxy = max(y)
h = barplot(y, las=2, col=colours, border=NA, ylim=c(0,maxy*1.5), space=.1, cex.names=0.6, names.arg=NULL, #names.arg=xstr, 
            ylab="Corr. mutation counts",main=paste("Nanoseq","\nobs_cos_sim=",round(obs_cosine_sim,2),"\nexp_cos_sim=",round(expected_cosine_sim,2)," [CI95: ",round(lci_95,2),"-",round(uci_95,2),"]",sep=""),cex.main=.7)
for (j in 1:length(sub_vec)) {
 	xpos = h[c((j-1)*16+1,j*16)]
    #rect(xpos[1]-0.5, maxy*1.2, xpos[2]+0.5, maxy*1.3, border=NA, col=colours[j*16])
 	#text(x=mean(xpos), y=maxy*1.3, pos=3, label=sub_vec[j])
}    
	

