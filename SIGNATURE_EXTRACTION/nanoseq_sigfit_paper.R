# NANOSEQ SIGNATURE ANALYSIS
# 2020/05/21

library(sigfit)
library(stringr)
#setwd("~/Desktop/nanoseq_sigfit/")


# A) SIGNATURE EXTRACTION
###########################

# Load counts and discard urothelium, bladder and colibactin-exposed samples
profiles = read.table("nanoseq.triprofiles.ok.tsv", sep="\t", header=T, as.is=T)
rownames(profiles) = profiles$SampleID
counts = round(profiles[!grepl("(PD37449)|(urothel)|(bladder)", profiles$SampleID), -(1:2)])
rownames(counts) = gsub(" nanoseq", "_nanoseq", gsub("nanoseqv", "nanoseq", rownames(counts)))


# Merge counts by tissue and sequencing type
types = c("cordblood_nanoseq", "cordblood_stdseq", "grans_nanoseq",
          "(adultblood_stdseq)|(colonies_stdseq)",
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
signatures = retrieve_pars(fit.merged.3, "signatures")$mean
#exposures = retrieve_pars(fit.merged.3, "exposures")$mean



# B) SIGNATURE EXPOSURE ANALYSES
##################################

# Load sample info table, reload counts
sample.data = read.table("Table_NanoSequenced_samples.txt", sep="\t", header=T, as.is=T, check.names=F)
profiles = read.table("nanoseq.triprofiles.ok.tsv", sep="\t", header=T, as.is=T)
rownames(profiles) = gsub(" nanoseq", "_nanoseq", gsub("nanoseqv", "nanoseq", profiles$SampleID))
profiles = profiles[, -(1:2)]


# Fit the extracted signatures to each sample
fit.sample.3 = fit_signatures(profiles, signatures, iter=5000, seed=0xC0FFEE, cores=4)

plot_all(fit.sample.3, out_path="SignatureFitting_AllSamples",
         exp_cex_names=0.9, exp_margin_bottom=14)


# Convert exposures to mutations per cell
# Burden per cell = num_muts / num_sites_called * genome_size
genome.size = 5722652910
neurons.counts = rowSums(profiles[grep("neurons", rownames(profiles)), ])
neurons.idx = match(str_split_fixed(names(neurons.counts), "_", 2)[, 1], sample.data$`Donor id`)
stopifnot(!any(is.na(neurons.idx)))

neurons.burden = neurons.counts / sample.data$`Number of sites called`[neurons.idx] * genome.size
neurons.expos = retrieve_pars(fit.sample.3, "exposures")$mean[names(neurons.counts), ]
stopifnot(identical(names(neurons.burden), rownames(neurons.expos)))

neurons.burden.sigs = neurons.burden * neurons.expos
stopifnot(max(abs(rowSums(neurons.burden.sigs) - neurons.burden)) < 1e-10)

# Output results
save(sample.data, profiles, signatures, neurons.expos, neurons.burden, neurons.burden.sigs,
     file="Neurons_Exposures_Burden.RData")

cat("\nDone\n")
