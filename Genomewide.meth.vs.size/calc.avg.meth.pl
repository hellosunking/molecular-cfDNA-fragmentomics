#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;
#use KSLIB::loadGenome qw/loadGenome/;
#use KSLIB::cigarUtil qw/fix_seq_from_CIGAR extract_size_from_CIGAR/;
#use KSLIB::Digitalize qw/digitalize/;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <in.log> [in.log]\n\n";
#	print STDERR "\nThis program is designed to \n\n";
	exit 2;
}

foreach my $file ( @ARGV ) {
	open IN, "$file" or die( "$!" );
	my ($C, $T) = (0, 0);
	while( <IN> ) {
		next if /^#/;
		chomp;
		my @l = split /\t/;	##chr10	424	170	78	157	56
		next unless $l[0] =~ /^chr\d+$/;
		$C += $l[2] + $l[4];
		$T += $l[3] + $l[5];
	}
	close IN;

	my $sid = $file;
	$sid = $1 if $sid =~ /\/(\d+).meth.log$/;
	print join("\t", $sid, $C, $T, $C/($C+$T)*100), "\n";
}

