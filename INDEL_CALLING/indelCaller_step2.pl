#!/usr/bin/perl -w

use strict;

# * Perl script
#     * Read all the indels, group by read bundle, mark sites in read bundles where we would believe an indel
#         * Add for each site if it is masked or not - when a final event overlaps a masked site, remove
#         * Reformat read bundle names to: RB:Z:1,11373,11754,GCT,TTA
#     * Extract read bundles from bam to create mini bam
#         * Remove duplicate flag
#     * samtools + bcftools to call indels
#     * read vcf output and decide which are ok and if no indel is called, warn!


my $botseq_bam_file = $ARGV[0]; # This is the nanoseq bam file
my $out_name        = $ARGV[1];
my $ref_genome      = "/lustre/scratch119/casm/team268im/fa8/data/hs37d5.fa"; # MODIFY accordingly

# Load indels:
my %indels;
my $counts;
while(<stdin>) { # input from "grep -vh BULK *bed | "   <==  those bed files are generated 
#                                                            by indelCaller.pl [step 1]
	chomp;
	my($chr,$pos0,$pos1,$info)             = (split(/\t/,$_  ))[0,1,2,3]    ;
	my($rb_id,$dp,$qpos,$context,$sw,$snp) = (split(/;/,$info))[0,1,2,3,4,5];
	my @tmp = split(/[:\-\|]/,$rb_id);
	#$rb_id  = "$tmp[0],".($tmp[1]+1).",$tmp[2],$tmp[3],$tmp[4]";
	$"      = ",";
	$rb_id  = "@tmp";
	$sw     =~ s/SW=//;
	$snp    =~ s/cSNP=//;
	$indels{$rb_id}->{$pos1}->{"dp"     }  = $dp     ;
	$indels{$rb_id}->{$pos1}->{"qpos"   }  = $qpos   ;
	$indels{$rb_id}->{$pos1}->{"context"}  = $context;
	$indels{$rb_id}->{$pos1}->{"sw"     }  = $sw     ;
	$indels{$rb_id}->{$pos1}->{"snp"    }  = $snp    ;
	$counts++;
}
print STDERR "$counts indel sites seen\n";
print STDERR scalar(keys(%indels)), " readbundles with indels\n";

open(O, ">$out_name.final.vcf") || die "noaslal\n";
my $header = "";
my $vcf_content = "";
foreach my $rb_id (keys %indels) {
	my($chr,$start,$end) = (split(/,/,$rb_id))[0,1,2];
	my(%good_sites,%bad_sites);
	foreach my $pos ( keys %{$indels{$rb_id}} ) {
		if($indels{$rb_id}->{$pos}->{"sw"} == 0 && $indels{$rb_id}->{$pos}->{"snp"} == 0) {
			$good_sites{$pos} = 1;
		} else {
			$bad_sites{$pos} = 0;
		}	
	}
	if(scalar(keys(%good_sites)) == 0) {
		next; # all indels overlap noisy or common SNP sites
	}
	$"=",";
	my @tmp = keys %{$indels{$rb_id}};
	my $positions = "@tmp";
	print STDERR "Step 1...\n";
	print STDERR "  samtools view -h $botseq_bam_file $chr:$start-$end | egrep \"($rb_id)|(^\@)\"  | perl /lustre/scratch119/casm/team268im/fa8/BOTSEQ/remove_dup_flag.pl | samtools view -bo $out_name.tmp1.bam -\n";
	`samtools view -h $botseq_bam_file $chr:$start-$end | egrep "($rb_id)|(^\@)"  | perl /lustre/scratch119/casm/team268im/fa8/BOTSEQ/remove_dup_flag.pl | samtools view -bo $out_name.tmp.bam -`;
	print STDERR "Step 2...\n";
	`samtools index $out_name.tmp.bam`;
	print STDERR "Step 3...\n";
	`samtools mpileup --no-BAQ  -d 250 -m 2 -F 0.5 -r $chr:$start-$end --BCF --output-tags DP,DV,DP4,SP -f $ref_genome -o $out_name.bcf $out_name.tmp.bam`;
	print STDERR "Step 4...\n";
	`bcftools index -f $out_name.bcf $out_name.indexed.bcf`;
	print STDERR "Step 5...\n";
	`bcftools call --skip-variants snps --multiallelic-caller --variants-only  -O v $out_name.bcf -o $out_name.tmp.vcf`;
	print STDERR "Step 6...\n";
	`bcftools norm -f $ref_genome $out_name.tmp.vcf > $out_name.tmp2.vcf`;
	#print "$rb_id:$positions:\n";
	my $get_header = 0;
	if($header eq "") {
		$get_header = 1;
	}
	print STDERR "\nChecking result\n";
	open(I, "$out_name.tmp2.vcf") || die "nolala\n";
	while(<I>) {
		if(/^#/) {
			if($get_header == 1) {
				print O $_;
				$header .= $_;
			}
			next;
		}
		print "RESULT: $_\n";
		my @tmp = split(/\t/,$_);
		$tmp[7] = "$tmp[7];RB=$rb_id";
		$tmp[6] = "PASS";
		# If there is overlap with SW or cSNP, flag it as MASKED
		for(my $i=$tmp[1]; $i<= ($tmp[1] + length($tmp[3])); $i++) {
			if(exists($indels{$rb_id}->{$i}) && ($indels{$rb_id}->{$i}->{"sw"} == 1 || $indels{$rb_id}->{$i}->{"snp"} == 1)) {
				$tmp[6] = "MASKED";
			}
		}
		$"="\t";
		#$vcf_content .= "@tmp\n";
		print O "@tmp";
	}
	close(I);
}
close O;

__END__

