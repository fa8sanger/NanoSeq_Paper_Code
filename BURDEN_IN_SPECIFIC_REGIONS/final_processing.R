
library(Biostrings)

args = commandArgs(trailingOnly=TRUE)
burden_regions = args[1]
out_name       = args[2]

genome_counts = vector()
genome_counts[c("ACA","ACC","ACG","ACT","ATA","ATC","ATG","ATT","CCA","CCC","CCG","CCT","CTA","CTC","CTG","CTT",
                "GCA","GCC","GCG","GCT","GTA","GTC","GTG","GTT","TCA","TCC","TCG","TCT","TTA","TTC","TTG","TTT")] = 
                c(115415924,66550070,14381094,92058521,117976329,76401029,105094288,142651503,105547494,75238490,
                  15801067,101628641,73791042,96335416,115950255,114180747,82414099,68090507,13621251,80004082,
                  64915540,54055728,86012414,83421918,112085858,88336615,12630597,126566213,119020255,112827451,
                  108406418,219915599);

results = read.table(burden_regions,row.names=1,header=T,sep="\t",stringsAsFactors=F)

colnames(results) = gsub("\\.",">",colnames(results))
subs = colnames(results)[35:130]
tri_prefix = sapply(subs,function(x) substr(x,1,3))
tris = names(genome_counts)

##########################################################################################
# Trinucleotide substitution profile
		
	tri_bg = as.vector(as.matrix(results[1,tris]))
	names(tri_bg) = tris
	tri_bg[setdiff(names(genome_counts),names(tri_bg))]=0
	
	colours        = rep(c("deepskyblue","black","firebrick2","gray","darkolivegreen3","rosybrown2"),each=16)
	sub_vec        = c("C>A","C>G","C>T","T>A","T>C","T>G")
	ctx_vec        = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
	full_vec       = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
	xstr           = paste(substr(full_vec,5,5), substr(full_vec,1,1), substr(full_vec,7,7), sep="")
	ordered_names  = paste(xstr,">",rep(c("A","G","T","A","C","G"),each=16),sep="")
	tmp_ = as.matrix(results[1,subs])
	tri_obs = as.vector(tmp_[1,])
	names(tri_obs) = colnames(tmp_)
	tri_obs[setdiff(ordered_names,names(tri_obs))]=0
	tri_obs        = tri_obs[ordered_names]
	
	pdf(width=9,height=7,file=paste(out_name,".trinuc-profile.pdf",sep=""))
	par(mfrow=c(2,1))
	# Difference to genome trinuc freqs
	difference_from_genome_wide <- as.vector((tri_bg/sum(as.numeric(tri_bg))) / (genome_counts/sum(genome_counts)))
	names(difference_from_genome_wide) = names(tri_bg)
	ratio2genome = difference_from_genome_wide
	barplot(difference_from_genome_wide,las=2,main="Obs/Genome trinuc. freqs",cex.names=.7,cex.lab=.7,col=rgb(0.2,0.4,.6,.5),ylim=c(0.0,max(difference_from_genome_wide)), xpd=FALSE)
	abline(h=1,col="red",lty=2,lwd=2)	
	tris = sapply(names(tri_obs),function(x) unlist(strsplit(x,">"))[1])
	y = tri_obs/difference_from_genome_wide[tris]; maxy = max(y)
	h = barplot(y, las=2, col=colours, border=NA, ylim=c(0,maxy*1.5), space=1, cex.names=0.6, names.arg=xstr, ylab="Corrected mutation counts",main="Corrected-to-genome mutation counts")
	for (j in 1:length(sub_vec)) {
	    xpos = h[c((j-1)*16+1,j*16)]
	    rect(xpos[1]-0.5, maxy*1.2, xpos[2]+0.5, maxy*1.3, border=NA, col=colours[j*16])
	    text(x=mean(xpos), y=maxy*1.3, pos=3, label=sub_vec[j])
	}    
	trint_subst_obs = tri_obs
	trint_onto_genome = y
	trint_bg = tri_bg
	 
	dev.off()
		
###################
# Correct and save mutation burden
	mut_burden = sum(trint_subst_obs)/sum(trint_bg)
	mut_burden_corrected = sum(trint_onto_genome)/sum(trint_bg)
	correction_factor = mut_burden_corrected/mut_burden
	
	burdens           = as.data.frame(matrix(nrow=2,ncol=5))
	rownames(burdens) = c("observed","corrected")
	colnames(burdens) = c("muts","total","burden","burden_lci","burden_uci")
	
	burdens["observed", "muts"  ]  = sum(trint_subst_obs)
	burdens["observed", "total" ]  = sum(trint_bg)
	burdens["observed", "burden"]  = sum(trint_subst_obs)/sum(trint_bg)
	burdens["corrected","muts"  ]  = sum(trint_onto_genome)
	burdens["corrected","total" ]  = sum(trint_bg)
	burdens["corrected","burden"]  = sum(trint_onto_genome)/sum(trint_bg)
	
	# Add confidence intervals:
	ci = poisson.test(sum(trint_subst_obs))$conf.int[1:2]
	#ci = binom.wilson(sum(trint_subst_obs),sum(trint_bg))[,c("lower","upper")]
	burdens["observed","burden_lci"] = ci[1]/sum(trint_bg)
	burdens["observed","burden_uci"] = ci[2]/sum(trint_bg)
	ci_corrected = ci*correction_factor
	#ci = binom.wilson(burdens["corrected","muts"],burdens["corrected","total"])[,c("lower","upper")]
	burdens["corrected","burden_lci"] = ci_corrected[1]/sum(trint_bg)
	burdens["corrected","burden_uci"] = ci_corrected[2]/sum(trint_bg)
	
	# Save to table:
	write.table(burdens,file=paste(out_name,".mut_burden.tsv",sep=""),sep="\t",quote=F,row.names=T,col.names=T)

