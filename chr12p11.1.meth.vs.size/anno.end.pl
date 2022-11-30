#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL
#

use strict;
use warnings;

if( $#ARGV < 1 ) {
	print STDERR "\nUsage: $0 <chr12p11.1.nucleosome.center.bed> <in.loci>\n\n";
#	print STDERR "\nThis program is designed to \n\n";
	exit 2;
}

my @nc;	##nucleosome center
open IN, "$ARGV[0]" or die( "$!" );
while( <IN> ) {
	chomp;
	my @l = split /\t/;	##chr12	34485766	34485767	Shendure
	push @nc, $l[2];
}
close IN;
push @nc, 0;
push @nc, 1e9;
@nc = sort {$a<=>$b} @nc;

my %dist;
open IN, "$ARGV[1]" or die( "$!" );
while( my $locus = <IN> ) {
	chomp( $locus );
	my $i;
	for($i=1; $i<=$#nc; ++$i) {
		last if $nc[$i] >= $locus;
	}
	if( abs($nc[$i] - $locus) < abs($nc[$i-1] - $locus) ) {
		print join("\t", $locus, $nc[$i],   $locus-$nc[$i]  ), "\n";
		$dist{$locus-$nc[$i]} ++;
	} else {
		print join("\t", $locus, $nc[$i-1], $locus-$nc[$i-1]), "\n";
		$dist{$locus-$nc[$i-1]} ++;
	}
}
close IN;

my @s = sort {$a<=>$b} keys %dist;
foreach my $i ( $s[0] .. $s[-1] ) {
	my $d = $dist{$i} || 0;
	print STDERR "$i\t$d\n";
}

