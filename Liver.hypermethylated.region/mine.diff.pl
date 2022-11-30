#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;

if( $#ARGV < 2 ) {
	print STDERR "\nUsage: $0 <T-cell.Msuite.CpG.meth.bedgraph> <Liver.Msuite.CpG.meth.bedgraph> <out.bed>\n\n";
	exit 2;
}

my $tcell = load_meth( $ARGV[0] );
my $liver = load_meth( $ARGV[1] );

open OUT, ">$ARGV[2]" or die( "$!" );

print STDERR "Mining CpGs\n";
foreach my $chr ( 1..22 ) {
	my $ref = $tcell->{ "chr$chr" };
	my $alt = $liver->{ "chr$chr" };

	foreach my $i ( keys %$ref ) {
		next unless exists $alt->{$i};

		if( $ref->{$i}<20 && $alt->{$i}>80 ) {
			print OUT join("\t", "chr$chr", $i-1, $i, $ref->{$i}, $alt->{$i}), "\n";
		}
	}
}
close OUT;

sub load_meth {
	my $file = shift;
	print STDERR "Loading $file\n";

	my %m;
	open IN, "$file" or die( "$!" );
	while( <IN> ) {
		chomp;
		my @l = split /\t/;	##chr1	10524	10525	100
		$m{$l[0]}->{$l[2]} = $l[3];
	}
	close IN;
	return \%m;
}

