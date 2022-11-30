#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;

if( $#ARGV < 2 ) {
	print STDERR "\nUsage: $0 <in.SE.sam> <out.file> <chr12.fa> [min.CpG=1] [max.Error=0]\n\n";
	exit 2;
}

my $minCpG = $ARGV[3] || 1;
my $maxERR = $ARGV[4] || 0;

## load sequence
my $fa = 'X';	## X is a placeholder
my $chr12fasta = $ARGV[2];
open IN, "$chr12fasta" or die $!;
<IN>;
while( <IN> ) {
	chomp;
	$fa .= uc $_;
}
close IN;

## load SAM file
open R, ">$ARGV[1]" or die( "$!" );

open IN, "$ARGV[0]" or die( "$!" );
while( my $r1 = <IN> ) {
	my @l = split /\t/, $r1, 12;
##SRR.513 0 chr12 34485223 0 101M * 0 0	AGGAGTGT AJJJJAJJF XG:Z:CT
	my $r2 = <IN>;
	my @k = split /\t/, $r2, 12;

	## this should be ensured by the previous program
	## which also filters out the unpaired reads
	if( $l[3] > $k[3] ) {
		print STDERR "Error pair: $l[0]!\n";
		next;
	}

	my ($fSeq1, $fQ1) = fix_seq_from_CIGAR($l[9], $l[10], $l[5]);
	my ($fSeq2, $fQ2) = fix_seq_from_CIGAR($k[9], $k[10], $k[5]);
#	print STDERR "$l[0]\n$l[3]\t$fSeq1\n$k[3]\t$fSeq2\n";

	## merge read 1 and read 2
	my $merge;
	if( $l[3]+length($fSeq1) < $k[3] ) {	## No overlap
		$merge = $fSeq1;
		my $Nholder = $k[3] - $l[3] - length($fSeq1);
		$merge .= 'N' x $Nholder;
		$merge .= $fSeq2;
	} else {
		$merge = substr( $fSeq1, 0, $k[3]-$l[3] );
		my $offset = $k[3]-$l[3];
		my $cmpBit = $l[3]+length($fSeq1) - $k[3];
		for( my $j=0; $j!=$cmpBit; ++$j ) {
			if( substr($fQ1, $j+$offset, 1) gt substr($fQ2, $j, 1) ) {
				$merge .= substr($fSeq1, $j+$offset, 1);
			} else {
				$merge .= substr($fSeq2, $j, 1);
			}
		}
		$merge .= substr( $fSeq2, $cmpBit );
	}
#	print STDERR "$merge\n\n";

	my ($CpG, $M, $U, $N) = (0, 0, 0, 0);
#	print STDERR "$l[0]\n$merge\n";
	if( $l[-1] =~ /XG:Z:CT/ ) {	## watson chain
		my $ref = substr( $fa, $l[3], length($merge)+1 );
		my $i = -1;
		while( ($i=index($ref, "CG", $i+1)) >= 0 ) {
			++ $CpG;
			my $s = substr( $merge, $i, 1 );
			if( $s eq 'C' ){ ++ $M; }
			elsif( $s eq 'T' ){ ++ $U; }
			else{ ++ $N; }
		}
	} elsif( $l[-1] =~ /XG:Z:GA/ ) { ## crick chain
		my $ref = substr( $fa, $l[3], length($merge)+1 );
		my $i = -1;
		while( ($i=index($ref, "CG", $i+1)) >= 0 ) {
			++ $CpG;
			my $s = substr( $merge, $i+1, 1 );
			if( $s eq 'G' ){ ++ $M; }
			elsif( $s eq 'A' ){ ++ $U; }
			else{ ++ $N; }
		}
	} else {	## unknown data structure
		print STDERR "Error record: $l[0]!\n";
		next;
	}

	next if $CpG < $minCpG || $N > $maxERR;
	## the output is NOT in BED format as both start/end are 1-based
	print R join("\t", $l[0], $l[3], $l[3]+length($merge)-1,
			$CpG, $M, $U, $N), "\n";
}
close IN;

close R;

sub fix_seq_from_CIGAR {
	my $seq   = shift;
	my $qual  = shift;
	my $cigar = shift;

	## discard hard-clip where the original sequence does not exist in $seq
	$cigar =~ s/^\d+H//;
	$cigar =~ s/\d+H$//;

	if( $cigar !~ /^\d+M$/ ) {	## there are indels here
		$cigar =~ s/([MID])/$1:/g;
		my @info = split /:/, $cigar;
		my ($mseq, $mqual) = ('', '');	## modified seq and qual
		my $curr = 0;
		foreach my $m ( @info ) {
			if( $m =~ /^(\d+)I$/ || $m =~ /^(\d+)S$/ ) {	## insertion/soft-clip: CAUTION!!! I WILL DISCARD IT!!!
				$curr += $1;
			} elsif ( $m =~ /^(\d+)D$/ )	{ ## deletion: CAUTION!!! I WILL ADD 'D' to indicate deletion !!!
				$mseq .= 'D' x $1;
				$mqual.= '!' x $1;
			} elsif( $m =~ /^(\d+)M$/ ) {
				$mseq .= substr($seq,  $curr, $1);
				$mqual.= substr($qual, $curr, $1);
				$curr += $1;
			}
		}
		return ($mseq, $mqual);
	} else {	## nothing to update
		return ($seq, $qual);
	}
}

