#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;

if( $#ARGV < 1 ) {
	print STDERR "\nUsage: $0 <in.meth.bed> <out.prefix>\n\n";
#	print STDERR "\nThis program is designed to \n\n";
	exit 2;
}

open H, ">$ARGV[1].high.bed" or die( "$!" );
open L, ">$ARGV[1].low.bed" or die( "$!" );

open IN, "$ARGV[0]" or die( "$!" );
while( <IN> ) {
	chomp;
	my @l = split /\t/;	##chr2	20051750	20051826	SRR8932917.sra.15	4	0	4	0
	if( $l[4]!=0 && $l[7]==0 ) {	##this read covers some CpGs, and there is NO errors or SNPs
		if( $l[5]==0 ) {
			print L "$_\n";
		} elsif( $l[6] == 0 ) {
			print H "$_\n";
		}
	}
}
close IN;

close H;
close L;

