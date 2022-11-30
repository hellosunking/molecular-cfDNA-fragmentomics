#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;
#use KSLIB::loadGenome qw/loadGenome/;
#use KSLIB::cigarUtil qw/fix_seq_from_CIGAR extract_size_from_CIGAR/;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <in.PE.sam>\n\n";
#	print STDERR "\nThis program is designed to \n\n";
	exit 2;
}

my %sam;
open IN, "$ARGV[0]" or die( "$!" );
while( <IN> ) {
	my @l = split /\t/, $_, 2;
	push @{$sam{$l[0]}}, $_;
}
close IN;

foreach my $sid ( sort keys %sam ) {
	my $records = $sam{$sid};
	next if $#$records != 1;	## read1 or 2 is NOT within the target region
	my @r1 = split /\t/, $records->[0];
	my @r2 = split /\t/, $records->[1];
	if( $r1[3] <= $r2[3] ) {
		print $records->[0], $records->[1];
	} else {
		print $records->[1], $records->[0];
	}
}

