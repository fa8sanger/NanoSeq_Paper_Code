library(sigfit)
library(stringr)
library(lme4)

library(deconstructSigs)
genome_size = sum(tri.counts.genome)*2

setwd("SIGNATURE_RATES")

##########################################################################################
# Extended Data Figure 7, panel e, Accumulation of mutations with age for each signature 
# in neurons

nanoseq = read.table("../DATA/rates.tsv",header=T,sep="\t",row.names=1)
nanoseq = nanoseq[which(nanoseq$cell_type=="neurons"),]
profiles = read.table("../DATA/profiles.tsv", sep="\t", header=T, row.names=1)
profiles = profiles[grep("neurons",rownames(profiles)),]

# Fit the extracted signatures to each sample
fit.sample.3 = fit_signatures(profiles, signatures, iter=5000, seed=0xC0FFEE, cores=4)

# Convert exposures to mutations per cell
# Burden per cell = num_muts / num_sites_called * genome_size
neurons.counts = rowSums(profiles[grep("neurons", rownames(profiles)), ])
neurons.burden = neurons.counts / nanoseq[names(neurons.counts),"total"] * genome_size
neurons.expos = retrieve_pars(fit.sample.3, "exposures")$mean[names(neurons.counts), ]
stopifnot(identical(names(neurons.burden), rownames(neurons.expos)))
neurons.burden.sigs = neurons.burden * neurons.expos
stopifnot(max(abs(rowSums(neurons.burden.sigs) - neurons.burden)) < 1e-10)

nanoseq$SigA = neurons.burden.sigs[,1]
nanoseq$SigB = neurons.burden.sigs[,2]
nanoseq$SigC = neurons.burden.sigs[,3]

par(mfrow=c(2,3))
modelSigA = lm(SigA ~ age, data=nanoseq)
modelSigA_disease = lm(SigA ~ age + disease_state, data=nanoseq)
summary(modelSigA)
summary(modelSigA_disease)
plot(nanoseq$age, nanoseq$SigA, pch=20,xlab="Age",ylab="SigA",
     xlim=c(0,105),ylim=c(0,max(nanoseq$SigA)),
     main="SigA (p(Age)=0.00833; p(AD)=0.01087)")
points(nanoseq[which(nanoseq$disease_state=="AD"),"age"],nanoseq[which(nanoseq$disease_state=="AD"),"SigA"],pch=20,col="red")
legend("topleft",legend=c("healthy","AD"),col=c("black","red"),pch=20)
abline(coef=coef(modelSigA),col="orange",lty=2)

modelSigB = lm(SigB ~ age, data=nanoseq)
modelSigB_disease = lm(SigB ~ age + disease_state, data=nanoseq)
summary(modelSigB)
summary(modelSigB_disease)
plot(nanoseq$age, nanoseq$SigB, pch=20,xlab="Age",ylab="SigB",
     xlim=c(0,105),ylim=c(0,max(nanoseq$SigB)),
     main="SigB (p(Age)=0.00457; p(AD)=0.03334)")
points(nanoseq[which(nanoseq$disease_state=="AD"),"age"],nanoseq[which(nanoseq$disease_state=="AD"),"SigB"],pch=20,col="red")
abline(coef=coef(modelSigB),col="orange",lty=2)

modelSigC = lm(SigC ~ age, data=nanoseq)
modelSigC_disease = lm(SigC ~ age + disease_state, data=nanoseq)
summary(modelSigC)
summary(modelSigC_disease)
plot(nanoseq$age, nanoseq$SigC, pch=20,xlab="Age",ylab="SigC",
     xlim=c(0,105),ylim=c(0,max(nanoseq$SigC)),
     main="SigC (p(Age)=1.8e-08; p(AD)=0.620)")
points(nanoseq[which(nanoseq$disease_state=="AD"),"age"],nanoseq[which(nanoseq$disease_state=="AD"),"SigC"],pch=20,col="red")
abline(coef=coef(modelSigC),col="orange",lty=2)

##########################################################################################
# Extended Data Figure 7, panel f, Accumulation of mutations with age for each signature
# in smooth muscle

nanoseq = read.table("../DATA/rates.tsv",header=T,sep="\t",row.names=1)
nanoseq = nanoseq[grep("muscle",nanoseq$cell_type),]
profiles = read.table("../DATA/profiles.tsv", sep="\t", header=T, row.names=1)
profiles = profiles[grep("muscle",rownames(profiles)),]

# Fit the extracted signatures to each sample
fit.sample.3.muscle = fit_signatures(profiles, signatures, iter=5000, seed=0xC0FFEE, cores=4)

# Convert exposures to mutations per cell
# Burden per cell = num_muts / num_sites_called * genome_size
muscle.counts = rowSums(profiles[grep("muscle", rownames(profiles)), ])
muscle.burden = muscle.counts / nanoseq[names(muscle.counts),"total"] * genome_size
muscle.expos = retrieve_pars(fit.sample.3.muscle, "exposures")$mean[names(muscle.counts), ]
stopifnot(identical(names(muscle.burden), rownames(muscle.expos)))
muscle.burden.sigs = muscle.burden * muscle.expos
stopifnot(max(abs(rowSums(muscle.burden.sigs) - muscle.burden)) < 1e-10)

nanoseq$SigA = muscle.burden.sigs[,1]
nanoseq$SigB = muscle.burden.sigs[,2]
nanoseq$SigC = muscle.burden.sigs[,3]

nanoseq$donor = substr(rownames(nanoseq),1,7)

modelSigA = lmer(SigA ~ age + (age -1|donor), data=nanoseq)
modelSigA_origin = lm(SigA ~ age + cell_type, data=nanoseq)
m1  = lmer(SigA ~ age + (age -1|donor), data=nanoseq)
m0  = lmer(SigA ~ 1   + (age -1|donor), data=nanoseq)
anova(m0,m1) # pvalue for age
m0  = lmer(SigA ~ age + (age -1|donor), data=nanoseq)
m1  = lmer(SigA ~ age + (age -1|donor) + (age -1|cell_type), data=nanoseq)
anova(m0,m1) # pvalue for cell_type
summary(modelSigA)
summary(modelSigA_origin)
plot(nanoseq$age, nanoseq$SigA, pch=20,xlab="Age",ylab="SigA",
     xlim=c(0,105),ylim=c(0,max(nanoseq$SigA)),
     main="SigA (p(Age)=0.01709; p(Origin)=0.1762)")
points(nanoseq[which(nanoseq$cell_type=="smuscleC"),"age"],nanoseq[which(nanoseq$cell_type=="smuscleC"),"SigA"],pch=20,col="red")
legend("topleft",legend=c("Bladder","Colon"),col=c("black","red"),pch=20)
abline(coef=fixef(modelSigA),col="orange",lty=2)

modelSigB = lmer(SigB ~ age + (age -1|donor), data=nanoseq)
modelSigB_origin = lm(SigB ~ age + cell_type, data=nanoseq)
m1  = lmer(SigB ~ age + (age -1|donor), data=nanoseq)
m0  = lmer(SigB ~ 1   + (age -1|donor), data=nanoseq)
anova(m0,m1) # pvalue for age
m0  = lmer(SigB ~ age + (age -1|donor), data=nanoseq)
m1  = lmer(SigB ~ age + (age -1|donor) + (age -1|cell_type), data=nanoseq)
anova(m0,m1) # pvalue for cell_type
summary(modelSigB)
summary(modelSigB_origin)
plot(nanoseq$age, nanoseq$SigB, pch=20,xlab="Age",ylab="SigB",
     xlim=c(0,105),ylim=c(0,max(nanoseq$SigB)),
     main="SigB (p(Age)=0.002108; p(Origin)=1)")
points(nanoseq[which(nanoseq$cell_type=="smuscleC"),"age"],nanoseq[which(nanoseq$cell_type=="smuscleC"),"SigB"],pch=20,col="red")
abline(coef=fixef(modelSigB),col="orange",lty=2)

modelSigC = lmer(SigC ~ age + (age -1|donor), data=nanoseq)
modelSigC_origin = lm(SigC ~ age + cell_type, data=nanoseq)
m1  = lmer(SigC ~ age + (age -1|donor), data=nanoseq)
m0  = lmer(SigC ~ 1   + (age -1|donor), data=nanoseq)
anova(m0,m1) # pvalue for age
m0  = lmer(SigC ~ age + (age -1|donor), data=nanoseq)
m1  = lmer(SigC ~ age + (age -1|donor) + (age -1|cell_type), data=nanoseq)
anova(m0,m1) # pvalue for cell_type
summary(modelSigC)
summary(modelSigC_origin)
plot(nanoseq$age, nanoseq$SigC, pch=20,xlab="Age",ylab="SigC",
     xlim=c(0,105),ylim=c(0,max(nanoseq$SigC)),
     main="SigC (p(Age)4.854e-05; p(Origin)=1)")
points(nanoseq[which(nanoseq$cell_type=="smuscleC"),"age"],nanoseq[which(nanoseq$cell_type=="smuscleC"),"SigC"],pch=20,col="red")
abline(coef=fixef(modelSigC),col="orange",lty=2)
