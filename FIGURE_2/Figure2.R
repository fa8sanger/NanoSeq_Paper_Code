##########################################################################################
# Panel b - Linear regression rates colonies / grans 
# And tests using regression models

setwd("FIGURE_2")

library(lme4)
library(deconstructSigs)
library(bootpredictlme4)		
genome_size = sum(tri.counts.genome)*2

######################
# Load NanoSeq data
# And correct mutation burdens in cord blood - we estimated 26.2% embryonic 
# mutations are missed due to matched normal filtering	
	
	nanoseq = read.table("../DATA/rates.tsv",header=T,sep="\t")
	nanoseq = nanoseq[which(nanoseq$cell_type %in% c("grans","cordblood")),]
	nanoseq = nanoseq[-which(nanoseq$sampleID %in% c("PD43976-Botseq","PD43976-Nanoseq")),]
	nanoseq = nanoseq[grep("botseq",nanoseq$sampleID,invert=T),]
	
	# Apply the correction for cord blood:
	nanoseq_correction = 1/(1-0.26) # 26.2% percentage of mutations lost. I need to multiply by nanoseq_correction
	nanoseq$mut_rate_final = nanoseq$burden
	nanoseq$lci_final      = nanoseq$burden_lci     
	nanoseq$uci_final      = nanoseq$burden_uci     
	nanoseq[which(nanoseq$age==0),"mut_rate_final"] =   nanoseq[which(nanoseq$age==0),"burden"    ] * nanoseq_correction
	nanoseq[which(nanoseq$age==0),"lci_final"     ] =   nanoseq[which(nanoseq$age==0),"burden_lci"] * nanoseq_correction
	nanoseq[which(nanoseq$age==0),"uci_final"     ] =   nanoseq[which(nanoseq$age==0),"burden_uci"] * nanoseq_correction
	
	# Translate to number of mutations per cell:
	nanoseq$muts_per_cell  = nanoseq$mut_rate_final* genome_size 
	nanoseq$lci_final      = nanoseq$lci_final     * genome_size 
	nanoseq$uci_final      = nanoseq$uci_final     * genome_size 
		
	# Remove cord blood from the lineal regression
	nanoseq_for_lm = nanoseq[which(nanoseq$age>0),]
	nanoseq_for_lm$donor = as.character(nanoseq_for_lm$age) # age works to identify different donors
	
	# Create donor column for mixed effect model:
	nanoseq_for_lm$donor = substr(nanoseq_for_lm$sampleID,1,7)
	
	# Lineal regression with random intercepts per donor
	nanoseq_lmer = lmer ( muts_per_cell ~ age + (age - 1|donor), data=nanoseq_for_lm, REML=F)
	nanoseq_coef = fixef(nanoseq_lmer) 
	confint(nanoseq_lmer)
	# Intercept 142.06 [CI95% -115.26, 414.20]
	# Slope 20.49 [CI95% 15.75,  25.13]
	
	# Calculate conf. intervals with bootstrapping (bootpredictlme4) - for plotting purposes
	ages = seq(from=0,to=100,by=7)
	nsim = 1000
	nanoseq_lmer_simul = sapply(ages,function(x) predict(nanoseq_lmer, newdata=data.frame(age=x), 
	                            re.form=NA, se.fit=TRUE, nsim=nsim)$ci.fit[,1])

##########################################################################################
# Load HSC/MPP colonies data
# These burdens have been corrected (see Methods)

	colonies = read.table("../DATA/colonies_rates.tsv",sep="\t",header=T)
	# Keep only samples with at least 200 Mbp of final coverage:
	colonies = colonies[which(colonies$c20_nanoseq_genome>200e6),]
	# Remove cord blood from regression:
	colonies = colonies[which(colonies$age>0),]
	# Translate into mutations per cell:
	colonies$muts_per_cell = colonies$corr_rate * genome_size
	colonies$donor = colonies$file

	# Regression:
	colonies_lmer = lmer(muts_per_cell ~ age + (age - 1|donor),data=colonies,REML=F)
	colonies_coef = fixef(colonies_lmer)
	confint(colonies_lmer)
	# Intercept 123.20 [CI95% 27.25, 218.61]
	# Slope 19.83 [CI95% 18.20, 21.45]

	# Calculate conf. intervals with bootstrapping (bootpredictlme4) - for plotting purposes
	ages = seq(from=0,to=100,by=7)
	nsim = 1000
	colonies_lmer_simul = sapply(ages,function(x) predict(colonies_lmer, newdata=data.frame(age=x), 
	                             re.form=NA, se.fit=TRUE, nsim=nsim)$ci.fit[,1])


##########################################################################################
# Test if the slopes for colonies and grans are significantly different:

	colonies$type = "colonies"
	nanoseq_for_lm$type = "grans"
	both_df = rbind(colonies[,c("age","muts_per_cell","donor","type")],
	                nanoseq_for_lm[,c("age","muts_per_cell","donor","type")])
	interaction_model = lmer(muts_per_cell ~ age*type + (age - 1|donor),data=both_df,REML=F)
	# age:typegrans: 0.33 [CI95% -4.58, 5.01]
	
	# Test if the interaction is significant (i.e. if the slopes are significantly different):
	nointeraction_model = lmer(muts_per_cell ~ age + type + (age - 1|donor),data=both_df,REML=F)
	anova(interaction_model,nointeraction_model)
	# p=0.92

##########################################################################################
# Regression combining grans and colonies:

	colonies$type = "colonies"
	nanoseq_for_lm$type = "grans"
	both_df = rbind(colonies[,c("age","muts_per_cell","donor","type")],
	                nanoseq_for_lm[,c("age","muts_per_cell","donor","type")])
	both_lmer = lmer(muts_per_cell ~ age + (age - 1|donor) + type, data=both_df,REML=F)
	fixef(both_lmer); confint(both_lmer)
	# Fixed effects:
	# Intercept 121.79 [CI95% 30.82, 211.69]
	# Slope (age) 19.86 [CI95% 18.34, 21.38]
	# typegrans 50.71 [CI95% -14.68, 119.59]
	
	both_lmer_notype = lmer(muts_per_cell ~ age + (age - 1|donor), data=both_df,REML=F)
	# Intercept 127.39 [CI95% 36.88, 217.20]
	# Slope (age) 19.83 [CI95% 18.30, 21.36]

	# Test if type (colonies/grans) significantly improves the model:
	anova(both_lmer,both_lmer_notype)
	# p=0.12
	

##########################################################################################
# And plot panel b:

	colonies = read.table("../DATA/colonies_rates.tsv",sep="\t",header=T)
	# Keep only samples with at least 200 Mbp of final coverage:
	colonies = colonies[which(colonies$c20_nanoseq_genome>200e6),]

		par(mgp=c(2, .7, 0))
		par(mar=c(3,3,3,3))
		plot(NULL,NULL,xlim=c(0,90),ylim=c(0,3000),ylab="Number of mutations per cell",xlab="Age",frame.plot=F)
		polygon(x=c(ages, rev(ages)), y=c(nanoseq_lmer_simul[1,],  rev(nanoseq_lmer_simul[2,])),  col=adjustcolor("firebrick",alpha.f=0.2), border=NA)
		polygon(x=c(ages, rev(ages)), y=c(colonies_lmer_simul[1,], rev(colonies_lmer_simul[2,])), col=adjustcolor("darkcyan",alpha.f=0.2), border=NA)

		means = vector()
		files = unique(colonies$file)
		all_results_emily = data.frame()
		ages = vector()
		ages[files]  = colonies[!duplicated(colonies$file),"age"]
		for(file in files) {
		    res = colonies[which(colonies$file == file),]
		    res$donor = file
		    res$age = ages[file]
		    cat(median(res$c20_genome),"\n")
			res$muts_per_cell = res$corr_rate * genome_size
			cat(file,":",ages[file],"\n")
		    boxplot(res$muts_per_cell,notch=T,
		            cex.lab=.7,cex.names=.7,cex.main=.7,cex.axis=.7,outcex=.4,
		            add=T,at=ages[file],frame.plot=F,axes=F,boxwex=4,col="darkcyan",
		            pch=20,outcol="gray")
		    means[file] = mean(res$muts_per_cell)
		    if(nrow(all_results_emily)==0) {
		    	all_results_emily = res
		    } else {
		    	all_results_emily = rbind(all_results_emily,res)
		    }
		}

		abline(fixef(colonies_lmer)[1], fixef(colonies_lmer)[2], col="darkcyan",lty=2)					
		segments(nanoseq$age,nanoseq$lci_final,nanoseq$age,nanoseq$uci_final)
		points(nanoseq$age,nanoseq$muts_per_cell,pch=20,cex=1.0,col="firebrick")
		coef = fixef(nanoseq_lmer)
		abline(coef[1], coef[2], col="firebrick",lty=2)

##########################################################################################
# Panel c - substitution profiles 59 year old donor
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


##########################################################################################
# Panel d - colonic crypts boxplots and comparison to LCM-Caveman crypts
	library(deconstructSigs)
	genome_size = sum(tri.counts.genome)*2

	rr = read.table("../DATA/rates_caveman.tsv",sep="\t",header=F)
	colnames(rr) = c("kk","file","sample","muts_caveman","c20_genome","muts","corr_muts","c20_nanoseq_genome","obs_rate","corr_rate")
	rr$muts_per_cell = rr$corr_rate * genome_size
	files = unique(rr$file)
	
	files = files[c(5,4,3)]
	ages = vector()
	ages[files] = c(36,38,54)
	
	samples = c("PD37449","PD34201","PD34200")

	# Now nanoseq:
	nanoseq = read.table("../DATA/rates.tsv",sep="\t",header=T)
	nanoseq = nanoseq[which(nanoseq$cell_type=="crypt"),]
	nanoseq$donor = substr(nanoseq$sampleID,1,7)
	nanoseq = nanoseq[which(nanoseq$donor %in% samples),]
	nanoseq$mut_rate_final = nanoseq$burden
	nanoseq$lci_final      = nanoseq$burden_lci     
	nanoseq$uci_final      = nanoseq$burden_uci     

	nanoseq$muts_per_cell  = nanoseq$mut_rate_final* genome_size 
	nanoseq$lci_final      = nanoseq$lci_final     * genome_size 
	nanoseq$uci_final      = nanoseq$uci_final     * genome_size 

	par(mgp=c(2, .7, 0))
	par(mar=c(3,3,3,3))
	means = vector()
	
	
	# and plot
	par(mar=c(8,3,3,3))
	box=boxplot(list("PD37449 (36F)"=c(rr[grep(samples[1],rr$file),"muts_per_cell"]),
	                 "PD34201 (38M)"=c(rr[grep(samples[2],rr$file),"muts_per_cell"]),
	                 "PD34200 (54M)"=c(rr[grep(samples[3],rr$file),"muts_per_cell"])),
	                 las=2,cex.names=.7,notch=T,ylab="Subs per cell",ylim=c(0,7000),main="Subs",
	                 col="darkcyan",pch=20,outcol="gray")
	nanoseq[grep(samples[1],nanoseq$tag),"x"] = 1
	nanoseq[grep(samples[2],nanoseq$tag),"x"] = 2
	nanoseq[grep(samples[3],nanoseq$tag),"x"] = 3
	points(rep(1,nrow(nanoseq[grep(samples[1],nanoseq$donor),])),nanoseq[grep(samples[1],nanoseq$donor),"muts_per_cell"],pch=20,col="firebrick",cex=1.3)
	points(rep(2,nrow(nanoseq[grep(samples[2],nanoseq$donor),])),nanoseq[grep(samples[2],nanoseq$donor),"muts_per_cell"],pch=20,col="firebrick",cex=1.3)
	points(rep(3,nrow(nanoseq[grep(samples[3],nanoseq$donor),])),nanoseq[grep(samples[3],nanoseq$donor),"muts_per_cell"],pch=20,col="firebrick",cex=1.3)




##########################################################################################
# Panel e
# Burden signatures 1 and 5
library(deconstructSigs)
profiles = read.table("../DATA/profiles.tsv",sep="\t",header=T,row.names=1)
profiles = profiles[grep("nanoseq",rownames(profiles)),]
profiles = profiles[grep("crypt",rownames(profiles)),]
colnames(profiles) = gsub("\\.",">",colnames(profiles))
exposures = as.data.frame(matrix(nrow=length(grep("colonic",rownames(profiles),value=T)),ncol=3))
colnames(exposures) = c("SBS1","SBS5","Colibactin") 
rownames(exposures) = grep("colonic",rownames(profiles),value=T)
for(s in grep("colonic",rownames(profiles),value=T)) {
	kk=whichSignatures(profiles,s, contexts.needed=T,
 	                  signatures.ref = "../DATA/one_five_colibactin.v3.tsv")
	cat(s,":","\n")
	print(kk$weights)
	exposures[s,] = kk$weights
}
par(mar=c(12,3,3,3))
barplot(t(as.matrix(exposures)),las=2,cex.names=.7,ylim=c(0,1.5),col=c("coral3","cornsilk4","darkgoldenrod"))
legend("top",legend=c("SBS1","SBS5","Colibactin") ,
       col=c("coral3","cornsilk4","darkgoldenrod"),pch=15,
       cex=.6,pt.cex=1)



# Plot regression signatures SBS1 + SBS5 (i.e. excluding colibactin)
genome_size = sum(tri.counts.genome)*2

nanoseq = read.table("../DATA/rates.tsv",header=T,sep="\t")
nanoseq = nanoseq[which(nanoseq$cell_type=="crypt"),]

nanoseq$mut_rate_final = nanoseq$burden
nanoseq$lci_final      = nanoseq$burden_lci     
nanoseq$uci_final      = nanoseq$burden_uci     

nanoseq$muts_per_cell  = nanoseq$mut_rate_final* genome_size 
nanoseq$lci_final      = nanoseq$lci_final     * genome_size 
nanoseq$uci_final      = nanoseq$uci_final     * genome_size 

nanoseq$sig1_muts = nanoseq$muts_per_cell * exposures[,"SBS1"]
nanoseq$sig5_muts = nanoseq$muts_per_cell * exposures[,"SBS5"]
nanoseq$sig1_mutslci = nanoseq$lci_final * exposures[,"SBS1"]
nanoseq$sig5_mutslci = nanoseq$lci_final * exposures[,"SBS5"]
nanoseq$sig1_mutsuci = nanoseq$uci_final * exposures[,"SBS1"]
nanoseq$sig5_mutsuci = nanoseq$uci_final * exposures[,"SBS5"]

nanoseq$donor = substr(nanoseq$sampleID,1,7)
nanoseq$sig1y5_muts    = nanoseq$muts_per_cell * (exposures[,"SBS1"]+exposures[,"SBS5"])
nanoseq$sig1y5_mutslci = nanoseq$lci_final     * (exposures[,"SBS1"]+exposures[,"SBS5"])
nanoseq$sig1y5_mutsuci = nanoseq$uci_final     * (exposures[,"SBS1"]+exposures[,"SBS5"])

# Regression and plot
	# With conf intervals and mixed-effects:
	library(lme4)
	model1 = lmer(sig1y5_muts ~ age + (age-1|donor),data=nanoseq)
	# Now without intercept! (because it was not significantly different from 0)
	model2 = lmer(sig1y5_muts ~ age + (age-1|donor) -1,data=nanoseq)
	plot(nanoseq$age,nanoseq$sig1y5_muts, type="n",
	     xlab="Age",
	     ylab="SBS1/SBS5 substitutions per cell",
	     ylim=c(0,max(nanoseq$sig1y5_mutsuci)),
	     main="Colonic crypts",
	     las=2,
	     xlim=c(0,80), frame.plot=F)
	library(bootpredictlme4)
	ages = seq(from=0,to=100,by=7)
	nsim = 1000
	rr = sapply(ages,function(x) predict(model1, newdata=data.frame(age=x), re.form=NA, se.fit=TRUE, nsim=nsim)$ci.fit[,1])		
	polygon(x=c(ages, rev(ages)), y=c(rr[1,], rev(rr[2,])), col="grey90", border=NA)
	segments(nanoseq$age,nanoseq$sig1y5_mutslci,nanoseq$age,nanoseq$sig1y5_mutsuci,col="black")
	points(nanoseq$age,nanoseq$sig1y5_muts,pch=20,col="firebrick")
	abline(coef=fixef(model1),col="firebrick",lty=2)



##########################################################################################
# Panel EDF6d: substitution profiles and cosine similarities:

profiles = read.table("../DATA/profiles.tsv",sep="\t",header=T,row.names=1)
nano_counts = list()
all_counts  = list()
all_sample_counts  = list() # separately by crypt
samples = c("PD37449","PD34201","PD34200")
for(sample in samples) {
	tmp = profiles[grep("nanoseq",rownames(profiles)),]
	tmp = tmp[grep("crypt",rownames(tmp)),]
	tmp = tmp[grep(sample,rownames(tmp)),]
	if(nrow(tmp) > 1) {
		tmp = apply(tmp,2,sum)
	} else {
		tmp = as.vector(as.matrix(tmp))
	}
	names(tmp) = gsub("\\.",">",colnames(profiles))
	nano_counts[[sample]] = tmp
	tmp = profiles[grep("stdseq",rownames(profiles)),]
	tmp = tmp[grep("crypt",rownames(tmp)),]
	tmp = tmp[grep(sample,rownames(tmp)),]
	if(nrow(tmp) > 1) {
		tmp = apply(tmp,2,sum)
	} else {
		tmp = as.vector(as.matrix(tmp))
	}
	names(tmp) = gsub("\\.",">",colnames(profiles))
	all_counts[[sample]] = tmp
	tmp = profiles[grep("stdseq",rownames(profiles)),]
	tmp = tmp[grep("crypt",rownames(tmp)),]
	for(i in c(1:nrow(tmp))) {
		j = as.vector(as.matrix(tmp[i,]))
		names(j) = gsub("\\.",">",colnames(profiles))
		all_sample_counts[[rownames(tmp)[i]]] = j
	}
}

# Now, for each sample, do the game of cosine similarites and simulations
for(sample in samples) {
		library(lsa)
		
		obss = exps = lcis = ucis = vector()
		
		obs_ = nano_counts[[sample]]
		ref_ = all_counts[[sample]]
		
		ref = obs = vector(length=96)
		names(ref) = names(obs) = names(ref_)
		
		obs = obs_[names(ref)]
		ref = ref_	
		
		# Calculate cosine similarity
		obs_cosine_sim = cosine(obs,ref)
		
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
		# Per crypt: a different way
		SIMUL_SIZE = 1000
		for(i in c(1:SIMUL_SIZE)) {
			# Simulate dataset of size round(sum(obs))
			size = round(sum(obs))
			ref_ = all_sample_counts[[sample(grep(sample,names(all_sample_counts),value=T),1)]]
			simul_ = sample(names(ref_),size=size,p=ref_/sum(ref_),replace=T)
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
		stdev_expected_cosine_sim = sd(cosines)  
		lci_95                    = sort(cosines)[round(SIMUL_SIZE*0.025)] 
		uci_95                    = sort(cosines)[round(SIMUL_SIZE*0.975)] 
		
		cat("Muts_qry\t",round(sum(obs)),"\n",sep="")
		cat("Muts_ref\t",round(sum(ref)),"\n",sep="")
		cat("obs_cos\t",obs_cosine_sim,"\n",sep="")
		cat("exp_cos\t",expected_cosine_sim,"\t+/-",stdev_expected_cosine_sim," [",lci_95,"-",uci_95,"]\n",sep="")
		
		colours        = rep(c("deepskyblue","black","firebrick2","gray","darkolivegreen3","rosybrown2"),each=16)
		sub_vec        = c("C>A","C>G","C>T","T>A","T>C","T>G")
		ctx_vec        = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
		full_vec       = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
		xstr           = paste(substr(full_vec,5,5), substr(full_vec,1,1), substr(full_vec,7,7), sep="")
		ordered_names  = paste(xstr,">",rep(c("A","G","T","A","C","G"),each=16),sep="")
		obs            = obs[ordered_names]
		ref            = ref[ordered_names]
		
		pdf(paste(sample,".profile+cos-sim.bycript.pdf",sep=""),width=10,height=6)
		par(mfrow=c(2,1))
		par(mar=c(3,6,3,3))
		par(mgp=c(3,1,0))
		y = obs
		maxy = max(y)
		h = barplot(y, las=2, col=colours, border=NA, ylim=c(0,maxy*1.5), space=1, cex.names=0.6, names.arg=xstr, 
		            ylab="Corrected mutation counts",main=paste("Nanoseq colonic crypts","\nobs_cos_sim=",round(obs_cosine_sim,2),"\nexp_cos_sim=",round(expected_cosine_sim,2)," [CI95: ",round(lci_95,2),"-",round(uci_95,2),"]",sep=""),cex.main=.7)
		for (j in 1:length(sub_vec)) {
		 	xpos = h[c((j-1)*16+1,j*16)]
		    rect(xpos[1]-0.5, maxy*1.2, xpos[2]+0.5, maxy*1.3, border=NA, col=colours[j*16])
		 	text(x=mean(xpos), y=maxy*1.3, pos=3, label=sub_vec[j])
		}    
		y = ref
		maxy = max(y)
		h = barplot(y, las=2, col=colours, border=NA, ylim=c(0,maxy*1.5), space=1, cex.names=0.6, names.arg=xstr, 
		            ylab="Corrected mutation counts",main="reference=Sigurgeir's calls")
		for (j in 1:length(sub_vec)) {
		 	xpos = h[c((j-1)*16+1,j*16)]
		    rect(xpos[1]-0.5, maxy*1.2, xpos[2]+0.5, maxy*1.3, border=NA, col=colours[j*16])
		 	text(x=mean(xpos), y=maxy*1.3, pos=3, label=sub_vec[j])
		}    
		dev.off()
}
