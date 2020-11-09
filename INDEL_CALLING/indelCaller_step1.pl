#!/usr/bin/perl -w

use strict;
$|=1;

if(scalar(@ARGV) < 3 || $ARGV[0] =~ /-h/) {
	print STDERR "Not all parameters have been correctly provided\n";
	print STDERR "Usage: perl table2bed_noErrors.pl min_fam_size filter_5_prime filter_3_prime < input_file.bed\n";
	print STDERR "       min_fam_size:   number, it can be 2 for 2+2, 3 for 3+3, etc\n";
	print STDERR "       filter_5_prime: number of bases to be trimmed from the 5', e.g. 10\n";
	print STDERR "       filter_3_prime: position of the trim from where to trim downstream bases.\n";
	print STDERR "                       For instance, if 10 bases are to be removed from the 3' of 150-bp reads,\n";
	print STDERR "                       then filter_3_prime should be specified as 141\n";
	exit;
}

my $MIN_SIZE_SUBFAM = $ARGV[0]; #Minimum number of family size 2 for 2+2, 3 for 3+3, etc
my $FILTER_5_PRIME  = $ARGV[1]; #Number of bases to be trimmed from 5'. Set as 0 if no filter is wanted. Example: 10 for the first 10 bases
my $FILTER_3_PRIME  = $ARGV[2]; #Number of bases to start trimming from. Set as 0 if no filter is wanted. For example: 140 for the 10 last bases of 150-bp reads
my $BULK_MIN_COV    = 16;

my %complement = ( "C" => "G", "G" => "C", "A" => "T", "T" => "A", "x" => "x", "?" => "?", "-" => "-", "N" => "N" ); 

while(<STDIN>) {
	next if(/^#/ || /chrom/); #in case it is receiving multiple csv files!
	chomp;

    my($chrom,$chromBeg,$chromEnd,$context,$commonSNP,$shearwater,$bulkASXS,$bulkNM,
       $bulkForwardA,$bulkForwardC,$bulkForwardG,$bulkForwardT,$bulkForwardIndel,
       $bulkReverseA,$bulkReverseC,$bulkReverseG,$bulkReverseT,$bulkReverseIndel,
       $dplxBreakpointBeg,$dplxBreakpointEnd,$dplxBarcode,$dplxOri,$dplxASXS,$dplxCLIP,$dplxNM,
       $dplxForwardA,$dplxForwardC,$dplxForwardG,$dplxForwardT,$dplxForwardIndel,
       $dplxReverseA,$dplxReverseC,$dplxReverseG,$dplxReverseT,$dplxReverseIndel,
       $dplxCQForwardA,$dplxCQForwardC,$dplxCQForwardG,$dplxCQForwardT,
       $dplxCQReverseA,$dplxCQReverseC,$dplxCQReverseG,$dplxCQReverseT,
       $bulkProperPair,$dplxProperPair) =split(/\t/,$_);
	   
	my $dplxForwardTotal = $dplxForwardA+$dplxForwardC+$dplxForwardG+$dplxForwardT+$dplxForwardIndel;
	my $dplxReverseTotal = $dplxReverseA+$dplxReverseC+$dplxReverseG+$dplxReverseT+$dplxReverseIndel;
	
	my $bulkForwardTotal = $bulkForwardA+$bulkForwardC+$bulkForwardG+$bulkForwardT+$bulkForwardIndel;
	my $bulkReverseTotal = $bulkReverseA+$bulkReverseC+$bulkReverseG+$bulkReverseT+$bulkReverseIndel;

	# Translate variable names (for compatibility with previous version):
    my $rbase  = (split(//,$context))[1];
    my $site   = $chromEnd;
    my $coord1 = $dplxBreakpointBeg;
    my $coord2 = $dplxBreakpointEnd;
    my $A      = $dplxForwardA;
    my $a      = $dplxReverseA;
    my $C      = $dplxForwardC;
    my $c      = $dplxReverseC;
    my $G      = $dplxForwardG;
    my $g      = $dplxReverseG;
    my $T      = $dplxForwardT;
    my $t      = $dplxReverseT;
	#my $r1     = $A+$C+$G+$T; # $fwdDP
	#my $r2     = $a+$c+$g+$t; # $revDP;
	my $r1     = $dplxForwardTotal;
	my $r2     = $dplxReverseTotal;
	my $dp     = ($r1+$r2)."[$r1,$r2]";
	my $chr    = $chrom;
	my $position_zerobased = $chromBeg;
	my $position_onebased  = $chromEnd;
	my $beg_breakpoint_coordinate = $dplxBreakpointBeg;
	my $end_breakpoint_coordinate = $dplxBreakpointEnd;

	my $bulk_fwd = $bulkForwardTotal;
	my $bulk_rev = $bulkReverseTotal;
	next if($bulk_fwd + $bulk_rev < $BULK_MIN_COV); # Bulk minimum coverage
	next if($dplxCLIP>0);
	next if($dplxNM > 20); # made very liberal to allow long indels. Check the impact!
	next if($dplxASXS < 50);

	if($r1 >= $MIN_SIZE_SUBFAM && $r2 >= $MIN_SIZE_SUBFAM) {
		my $bulktotal = $bulkForwardTotal+$bulkReverseTotal;
		if($dplxForwardIndel+$dplxReverseIndel < 0.9*($dplxForwardTotal+$dplxReverseTotal)) {
			next; # only interested in indels
		}
		my $qpos;
		my $orientation_type;
		my $qposF = $chromBeg-$dplxBreakpointBeg + 1;
		my $qposR = $dplxBreakpointEnd-$chromEnd;
		if($qposR < $qposF) {
			$qpos = $qposR;
			$orientation_type = "Rev";
		} else {
			$qpos = $qposF;
			$orientation_type = "Fwd";
		}
		if($FILTER_5_PRIME > 0) {
			if($qpos < $FILTER_5_PRIME) {
				next;
			}
		}
		if($FILTER_3_PRIME > 0) {
			if($qpos > $FILTER_3_PRIME) {
				next;
			}
		}
		my $dp = $r1+$r2;
		my $site_tags = "";
		$site_tags .= "$chr:$coord1-$coord2:$dplxBarcode;DP=$dp;QPOS=$qpos;";
		# If seen in the bulk, flag it:
		if($bulkForwardIndel+$bulkReverseIndel > 1) {
			$site_tags .= "BULK_SEEN($dplxForwardIndel+$dplxReverseIndel/$bulktotal);";
		}
		
		# Annotate the sequence context:
		my $signature;
		my $signature_trinuc;
		$signature_trinuc = $context    ;
		$signature_trinuc =~ s/\./$rbase/;
		$signature_trinuc .= ">indel";
		if($rbase =~ /[AGag]/) {
			$signature_trinuc = &reverse_signature($signature_trinuc);
		}
		$site_tags .= "$signature_trinuc;";
		print "$chr\t",$site-1,"\t$site\t$site_tags","SW=$shearwater;cSNP=$commonSNP;\n"; #BED Format incluyendo los site tags!!!!!!!!!!!!!
	}
}

sub reverse_signature {
	my $signature = uc($_[0]);
	my @tmp = split(//,$signature);
	$signature = $complement{$tmp[2]}.$complement{$tmp[1]}.$complement{$tmp[0]}.">"."indel";
	return $signature;
}




__END__


