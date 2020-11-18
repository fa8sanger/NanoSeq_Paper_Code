##########################################################################################
# Script to compare HSC/MPP colonies and differentiated granulocytes 
##########################################################################################

library(lme4)
library(deconstructSigs)
library(bootpredictlme4)		
genome_size = sum(tri.counts.genome)*2

setwd("BLOOD_REGRESSIONS")

##########################################################################################
# Load NanoSeq data
# And correct mutation burdens in cord blood - we estimated 26.3% embryonic 
# mutations are missed due to matched normal filtering	
	
	nanoseq = read.table("grans+cord.tsv",header=T,sep="\t")
	
	# Apply the correction for cord blood:
	nanoseq_correction = 1/(1-0.26) # 26.2% percentage of mutations lost. I need to multiply by nanoseq_correction
	nanoseq$mut_rate_final = nanoseq$mut_rate
	nanoseq$lci_final      = nanoseq$lci     
	nanoseq$uci_final      = nanoseq$uci     
	nanoseq[which(nanoseq$age==0),"mut_rate_final"] =   nanoseq[which(nanoseq$age==0),"mut_rate"] * nanoseq_correction
	nanoseq[which(nanoseq$age==0),"lci_final"     ] =   nanoseq[which(nanoseq$age==0),"lci"] * nanoseq_correction
	nanoseq[which(nanoseq$age==0),"uci_final"     ] =   nanoseq[which(nanoseq$age==0),"uci"] * nanoseq_correction
	
	# Translate to number of mutations per cell:
	nanoseq$muts_per_cell  = nanoseq$mut_rate_final* genome_size 
	nanoseq$lci_final      = nanoseq$lci_final     * genome_size 
	nanoseq$uci_final      = nanoseq$uci_final     * genome_size 
	
	# Remove non "NanoSeq" samples (standard BotSeqS and NanoSeq without ddBTPs):
	nanoseq = nanoseq[which(nanoseq$tag!="Botseq"&nanoseq$tag!="Nanoseq"),]
	
	# Remove cord blood from the lineal regression
	nanoseq_for_lm = nanoseq[grep("EM_",nanoseq$type,invert=T),]
	nanoseq_for_lm$donor = as.character(nanoseq_for_lm$age) # age works to identify different donors
	
	# Lineal regression with random intercepts per donor
	nanoseq_lmer = lmer ( muts_per_cell ~ age + (age - 1|donor), data=nanoseq_for_lm, REML=F)
	nanoseq_coef = fixef(nanoseq_lmer) 
	confint(nanoseq_lmer)
	# Intercept 140.46 [CI95% -117.56, 413.18]
	# Slope 20.54 [CI95% 15.79, 25.19]
	
	# Calculate conf. intervals with bootstrapping (bootpredictlme4) - for plotting purposes
	# ages = seq(from=0,to=100,by=7)
	# nsim = 1000
	# nanoseq_lmer_simul = sapply(ages,function(x) predict(mmm, newdata=data.frame(age=x), 
	#                             re.form=NA, se.fit=TRUE, nsim=nsim)$ci.fit[,1])



##########################################################################################
# Load HSC/MPP colonies data
# These burdens have been corrected (see Methods)

	colonies = read.table("colonies_rates.tsv",sep="\t",header=T)
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
	# ages = seq(from=0,to=100,by=7)
	# nsim = 1000
	# colonies_lmer_simul = sapply(ages,function(x) predict(mmmm, newdata=data.frame(age=x), 
	#                              re.form=NA, se.fit=TRUE, nsim=nsim)$ci.fit[,1])

	
##########################################################################################
# Test if the slopes for colonies and grans are significantly different:

	colonies$type = "colonies"
	nanoseq_for_lm$type = "grans"
	both_df = rbind(colonies[,c("age","muts_per_cell","donor","type")],
	                nanoseq_for_lm[,c("age","muts_per_cell","donor","type")])

	interaction_model = lmer(muts_per_cell ~ age*type + (age - 1|donor),data=both_df,REML=F)
	# age:typegrans: 0.30 [CI95% -4.50, 5.09]
	
	# Test if the interaction is significant (i.e. if the slopes are significantly different):
	nointeraction_model = lmer(muts_per_cell ~ age + type + (age - 1|donor),data=both_df,REML=F)
	anova(interaction_model,nointeraction_model)
	# p=0.90


##########################################################################################
# Regression combining grans and colonies:

	colonies$type = "colonies"
	nanoseq_for_lm$type = "grans"
	both_df = rbind(colonies[,c("age","muts_per_cell","donor","type")],
	                nanoseq_for_lm[,c("age","muts_per_cell","donor","type")])
	both_lmer = lmer(muts_per_cell ~ age + (age - 1|donor) + type, data=both_df,REML=F)
	fixef(both_lmer); confint(both_lmer)
	# Fixed effects:
	# Intercept 121.26 [CI95% 30.31, 211.14]
	# Slope (age) 19.86 [CI95% 18.35, 21.38]
	# typegrans 52.24 [CI95% -13.13, 121.05]
	
	both_lmer_notype = lmer(muts_per_cell ~ age + (age - 1|donor), data=both_df,REML=F)
	# Intercept 127.03 [CI95% 36.88, 217.20]
	# Slope (age) 19.84 [CI95% 18.31, 21.37]

	# Test if type (colonies/grans) significantly improves the model:
	anova(both_lmer,both_lmer_notype)
	# p=0.12
	

##########################################################################################

