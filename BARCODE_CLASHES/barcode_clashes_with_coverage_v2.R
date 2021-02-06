# Simulation of barcode clashes in nanoseq data as a function of the barcode frequency distribution

# 1. Multinomial simulation (it assumes exact same coverage at each site)

coverage = 10
barcodes = 500
aux = apply(rmultinom(n=1e4, size=coverage, prob=rep(1/barcodes,barcodes)), 2, function(x) c(sum(x==1), sum(x>=1)) )
sum(aux[1,])/sum(aux[2,]) # Estimate from the multinomial simulation
dbinom(x=1, size=coverage, p=1/barcodes) / pbinom(q=0, size=coverage, p=1/barcodes, lower.tail=F) # Exact binomial estimate

# 2. Poisson model (it assumes that molecules are randomly sampled, and so coverage has Poisson variation across sites)
#    This is more realistic than a multinomial model because mean genome coverage is not an integer value and can be <1.
#    A negative binomial model could be used learning the overdispersion (variation in coverage across sites) from nanoseq data.

coverage = 10
barcodes = 5
lambda = coverage/barcodes # Coverage per barcode

aux = sapply(1:1e4, function(j) { x = rpois(n = barcodes, lambda = lambda); c(sum(x==1), sum(x>=1)) })
sum(aux[1,])/sum(aux[2,]) # Estimate from the Poisson simulation
dpois(x=1, lambda=lambda) / ppois(q=0, lambda=lambda, lower.tail=F)

# 3. Poisson model with an input barcode frequency vector

bfreqs = read.table("200326_GlobalBarcodeFrequencies.csv", header=1, sep=",", stringsAsFactors=F)
f = rowMeans(as.matrix(bfreqs[,2:4]))
f = f / sum(f) # Average barcode frequencies across 3 datasets

coverage = 10 # Total coverage per site (number of molecules sampled)
lambda = coverage * f # Coverage per barcode (f = barcode frequency vector, summing to 1)

aux = sapply(1:1e4, function(j) { x = rpois(n = length(lambda), lambda = lambda); c(sum(x==1), sum(x>=1)) })
sum(aux[1,])/sum(aux[2,]) # Estimate from the Poisson simulation
# The weighting factor (how much each barcode weights in the mean) is proportional 
# to their relative coverage considering clashes (hence the use of the probability for families >=1)
relfootprint = ppois(q=0, lambda=lambda, lower.tail=F) 
sum(dpois(x=1, lambda=lambda) * relfootprint / ppois(q=0, lambda=lambda, lower.tail=F)) / sum(relfootprint)

# 4. Plotting the impact of coverage (assuming uniform -Poisson- coverage across sites)

cov_vec = 10^seq(-1,3,length.out = 100)
burden_subest = burden_subest_equibarcode = rep(NA, length(cov_vec))

for (j in 1:length(cov_vec)) {
    # Empirical barcode frequencies
    lambda = cov_vec[j] * f
    burden_subest[j] = sum(dpois(x=1, lambda=lambda) * relfootprint / ppois(q=0, lambda=lambda, lower.tail=F)) / sum(relfootprint)
    # Equal frequencies of all barcodes
    lambda = cov_vec[j] / 4096
    burden_subest_equibarcode[j] = sum(dpois(x=1, lambda=lambda) * relfootprint / ppois(q=0, lambda=lambda, lower.tail=F)) / sum(relfootprint)
}

dev.new(width = 5, height = 5)
plot(cov_vec, burden_subest, log="x", las=1, type="l", ylab="Burden subestimation", xlab="Molecules sampled per site (~coverage)")
lines(cov_vec, burden_subest_equibarcode, col="cadetblue")
abline(h = 0.99, col="grey")
dev.copy(pdf, file="barcode_clashes_burden_effect.pdf", width=5, height=5, useDingbats=F); dev.off(); dev.off()

