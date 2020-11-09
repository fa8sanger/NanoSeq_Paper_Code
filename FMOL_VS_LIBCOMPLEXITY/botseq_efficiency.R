
###########################################################################################
## Simulations to estimate the optimal bottleneck size given an estimate of 
## library complexity.
## Equations below are assuming 1 fmol sequenced equates to 1e8 molecules
## The relationship between fmol and molecules was derived empirically (see next section)
###########################################################################################

	# Truncated Poisson (only defined for x>=1): this models the duplicate family size distribution
	dpoistrunc = function(x, lambda) {
	    d = dpois(x, lambda) / (1-dpois(0, lambda))
	    d[x==0] = NA
	    return(d)
	}
	
	# The mean sequencing cost of a family is the mean of the zero-truncated Poisson: 
	#  lambda/(1-exp(-lambda))
	
	# Efficiency: The fraction of families of size>=2 divided by the average size (i.e. 
	#  sequencing cost) of a family, squared (as we need coverage from both strands)
	
	ratio = seq(0,20,by=0.1) # Sequencing output / library complexity
	
	efficiency = (ppois(q=2-0.1, lambda=ratio/2, lower.tail=F)/(1-dpois(0, ratio/2)))^2 / (ratio/(1-exp(-ratio)))
	
	# Duprate can be simply calculated using the mean of the family size (from the zero-truncated 
	#  Poisson): 1-1/mean = 1-(1-exp(-r))/r
	duprate = 1-(1-exp(-ratio))/ratio
	
	dev.new(width=9, height=3)
	par(mfrow=c(1,3))
	
	plot(ratio, duprate, ylim=c(0,1), type="l", xlab="Seq ratio (sequencing yield / lib complexity)", ylab="Duplicate rate", las=1)
	abline(h=1, lty=2)
	#lines(ratio, pmax(0, 1-1/ratio), col="palegreen3", lty=2) # Simplistic expectation (ratio=10 would tend to yield ~90% duplicate rates)
	
	plot(duprate, efficiency, type="l", xlim=c(0,1), xlab="Duplicate rate", ylab="Efficiency (duplex bp / total bp)", las=1)
	opt_duprate = duprate[which.max(efficiency)] # Optimal duplicate rate appears to be 0.805
	ind = which(efficiency > max(efficiency, na.rm=T)*0.8)
	abline(v=opt_duprate, col="red")
	abline(v=duprate[c(min(ind),max(ind))], col="red", lty=2)
	semiopt_duprate = duprate[c(min(ind),max(ind))]
	
	plot(ratio, efficiency, type="l", xlab="Seq ratio (sequencing yield / lib complexity)", ylab="Efficiency (duplex bp / total bp)", las=1)
	opt_ratio = ratio[which.max(efficiency)] # Optimal duplicate rate appears to be 5.1
	ind = which(efficiency > max(efficiency, na.rm=T)*0.8)
	abline(v=opt_ratio, col="red")
	abline(v=ratio[c(min(ind),max(ind))], col="red", lty=2)
	semiopt_ratio = ratio[c(min(ind),max(ind))]
	
	dev.copy(pdf, file="Nanoseq_efficiency_plots.pdf", width=9, height=3, useDingbats=F); dev.off(); dev.off()



###########################################################################################
## Linear relationship between library complexity and qPCR quantification
## Our empirical results indicate that 1 fmol ~ 1e8 molecules 
## (our size selection aims at 250-500 bps)
###########################################################################################

	## MySeq library (Aug 2019)

	qpcr = read.table("miseq_concentrations_Pete_Ellis_20190809.txt", header=1, sep="\t", stringsAsFactors=F)
	plot(qpcr$fmol, qpcr$unique_fragments, xlab="Library yield (fmol)", ylab="Unique fragments", las=1, pch=20, col="cadetblue4")
	model = lm(qpcr$unique_fragments ~ qpcr$fmol -1)
	abline(model, col="orange", lty=2)
	summary(model)
	coefficients(model) # slope (unique fragments per fmol): (9.8e7 unique fragments / fmol)
	confint(model) # Confidence intervals for the slope (CI95%: 9.5e7, 10.1e7)
	
	#################################################################
	# Using the median of the fragment/fmol relationship
	median(qpcr$unique_fragments/qpcr$fmol) # 10.27e7 instead of 9.8e7
	
	#################################################################
	# Based on the median and the regression we establish that the number of fragments 
	# per fmol is ~ 1e8 (with our size selection for fragments between 250-500 bps)
	#################################################################
	
	#################################################################
	# Calculating the volume required to do a deep NanoSeq experiment (equivalent to a 15x standard)
	# The 5.1 factor is the duplicate ratio estimated in the analytical simulation 
	qpcr$fmol_per_ul = qpcr$fmol / 12.5   # 12.5 ul taken
	coverage         = 15                 # coverage desired in genome fold units (e.g. 15x)
	fragmperfmol     = 1e8             
	qpcr$vol = round(coverage * 1e7 / (qpcr$fmol_per_ul * fragmperfmol * 5.1))	
	qpcr$vol_fromMiSeq = round(coverage * 1e7 * 12.5 / (qpcr$unique_fragments * 5.1), digits = 1)
	





