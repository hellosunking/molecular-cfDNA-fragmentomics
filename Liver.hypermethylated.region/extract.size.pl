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
	print STDERR "\nUsage: $0 <in.alleles> <out.prefix>\n\n";
#	print STDERR "\nThis program is designed to \n\n";
	exit 2;
}

open IN, "$ARGV[0]" or die( "$!" );
my %size;
my (%m, %u);
while( <IN> ) {
	chomp;
	my @l = split /\t/;	##chr1	1034052	1034193	CG	JF	CRR053121.sra.4676639	1034163	CT

	if( $l[-1] eq 'CT' ) {
		if( $l[3] eq 'CG' ) {
			$m{$l[5]} = 1;
			$size{$l[5]} = $l[2] - $l[1];
		} elsif( $l[3] eq 'TG' ) {
			$u{$l[5]} = 1;
			$size{$l[5]} = $l[2] - $l[1];
		} else {
			## discard
		}
	} else {
		if( $l[3] eq 'CG' ) {
			$m{$l[5]} = 1;
			$size{$l[5]} = $l[2] - $l[1];
		} elsif( $l[3] eq 'CA' ) {
			$u{$l[5]} = 1;
			$size{$l[5]} = $l[2] - $l[1];
		} else {
			## discard
		}
	}
}
close IN;

my (%sm, %su);
my ($m_all, $u_all) = ( 0, 0 );
foreach my $sid ( keys %m ) {
	next if exists $u{$sid};	## conflict result
	$sm{ $size{$sid} } ++;
	++ $m_all;
}
foreach my $sid ( keys %u ) {
	next if exists $m{$sid};	## conflict result
	$su{ $size{$sid} } ++;
	++ $u_all;
}

open OUT, ">$ARGV[1].methed.size" or die("$!");
my $cumu = 0;
foreach my $i ( 1..1000 ) {
	my $here = $sm{$i} || 0;
	$cumu += $here;
	print OUT join("\t", $i, $here, $cumu/$m_all), "\n";
}
close OUT;

open OUT, ">$ARGV[1].unmeth.size" or die("$!");
$cumu = 0;
foreach my $i ( 1..1000 ) {
	my $here = $su{$i} || 0;
	$cumu += $here;
	print OUT join("\t", $i, $here, $cumu/$u_all), "\n";
}
close OUT;

