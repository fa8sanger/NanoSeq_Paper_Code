#!/usr/bin/perl -w

##########################################################################################
# Fa8, 2020
# Script to identify which regions in a bam file have at least a given coverage 
# (20X in the manuscript). Reads with poor mapping quality (<10) and duplicates are 
# ignored.
#
# Requirements: samtools and bedtools
# 
# The script expects a bam file and a coverage threshold as arguments
# It also expects to find a file named as the bam file + "caveman_c.no.analysis.bed.gz", 
# which contains regions not analysed by Caveman
##########################################################################################

use strict;

my $bam = $ARGV[0];
my $cov = $ARGV[1]; # e.g. 20
my($sample_id);
if($bam =~ /\//) { # full path
	$sample_id = (split(/[\/\.]/,$bam))[1];
} else {
	$sample_id = (split(/[\/\.]/,$bam))[0];
}
print STDERR "bam=$bam; sample_id=$sample_id;\n";
my $no_analysis = "$sample_id.caveman_c.no.analysis.bed.gz";
my $human_genome = "human_genome.bed";
`zcat $no_analysis | subtractBed -a $human_genome -b stdin > $sample_id.analysed_genome.bed`;
open(O, ">$sample_id.min$cov.bed.tmp") || die "anñdkja\n";
open(I, "samtools view -h -q 10 -F 1024  $bam | genomeCoverageBed -dz -ibam stdin | ") || die "noalalla\n";
my($prev_chr,$prev_start,$first_chr,$first_start) = ("","","","");
while(<I>) {
	chomp;
	my($chr,$start,$depth) = (split(/\t/,$_))[0,1,2];
	if($prev_chr eq "") {
		if($depth >= $cov) {
			$first_chr   = $chr;
			$first_start = $start;
		}
		$prev_chr   = $chr;
		$prev_start = $start;
		next;
	} 
	if($depth >= $cov && $prev_chr eq $chr) {
		if($first_start eq "") {
			$first_chr   = $chr;
			$first_start = $start;
		}
	} else {
		if($depth < $cov && $first_start ne "") {
			print O "$chr\t$first_start\t",$prev_start+1,"\n";
			$first_chr   = "";
			$first_start = "";
		} elsif($prev_chr ne $chr) {
			$first_chr   = $chr;
			$first_start = $start;
		}
	}
	$prev_chr   = $chr;
	$prev_start = $start;
}
close(I);
close(O);

`intersectBed -a $sample_id.min$cov.bed.tmp -b $sample_id.analysed_genome.bed  > $sample_id.min$cov.bed`;
`/software/R-3.6.1/bin/Rscript /lustre/scratch119/casm/team268im/fa8/BOTSEQ/SCRIPTS/trinuc_counts_by_chunks.R $sample_id.min$cov.bed`;

__END__

