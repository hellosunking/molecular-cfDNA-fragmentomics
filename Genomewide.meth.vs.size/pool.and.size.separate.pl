#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <out.dir> <in.nameSrt.bam> [in.nameSrt.bam ...]\n\n";
	exit 2;
}

my $outDIR = shift;
my %cnt;
foreach my $file ( @ARGV ) {
	print STDERR "Loading $file\n";

	my %sam;
	open IN, "samtools view -@ 4 $file |" or die( "$!" );
	while( my $r1 = <IN> ) {
		my $r2 = <IN>;
		my @l = split /\t/, $r1, 10;
		##CRR053084.sra.11580582	163	chr1	10006	0	64M1I35M	=	10104	197	CTAACCCTA	AAFFFJJJJJJJJJJJJJ
		my $size = abs($l[8]);
		if( $size > 0 ) {
			$sam{$size} .= "$r1$r2";
		} else {
			print STDERR "ERROR: unknown size on line $.!\n";
		}
		$cnt{$size} ++;
	}
	close IN;

	print STDERR "Writing output\n";
	foreach my $s ( 50..500 ) {
		open OUT, ">>$outDIR/$s.sam" or die( "$!" );
		print OUT $sam{$s};
		close OUT;
	}
}

open OUT, ">$outDIR/size.stat" or die( "$!" );
foreach my $s ( sort {$a<=>$b} keys %cnt ) {
	print OUT "$s\t$cnt{$s}\n";
}

