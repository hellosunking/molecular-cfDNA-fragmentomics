#!/usr/bin/perl
#
# Author: Kun Sun @SZBL (sunkun@szbl.ac.cn)
#

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <pe.bed> [autosomeOnly=n] [fragLimit=inf]\n";
	print STDERR "\nThis program is designed to extract the fragment size from paired-end soap or sam/bam file.";
	print STDERR "\nBy setting the fragLimit, only the first fragLimit lines will be loaded (with a higher speed).\n\n";
	exit 2;
}

my $MAXVALUE = 2000;

my $file = $ARGV[0];

my $autosomeOnly = 0;
if( defined $ARGV[1] ) {
	if( $ARGV[1]=~/^Y(ES)?$/i ) {
		$autosomeOnly = 1;
		print STDERR "INFO: Autosome Only Mode is ON.\n";
	}
}

my $fragLimit = 'inf';
if( defined $ARGV[2] ) {
	if( $ARGV[2] > 0 ) {
		print STDERR "INFO: fragment limit is set to $ARGV[2].\n";
		$fragLimit = $ARGV[2] * 2;
	} elsif( $ARGV[2] == 0 ) {
		## do nothing, 0 stands for infinity
	} else {
		print STDERR "WARNING: Invalid fragment limit! IGNORED.\n";
	}
}

my @fraglen;
if( $ARGV[0] =~ /\.gz$/ ) {
	open IN, "zcat $ARGV[0] |" or die("$!");
} elsif( $ARGV[0] =~ /\.bz2$/ ) {
	open IN, "bzcat $ARGV[0] |" or die("$!");
} else {
	open IN, "$ARGV[0]" or die("$!");
}
my $all = 0;
while( <IN> ) {
	chomp;
	my @l = split /\t/;	##chr start end extra
	last if $. >= $fragLimit;
	next if $autosomeOnly && $l[0]!~/\d$/;
	my $len = $l[2] - $l[1];
#	if( $len < 0 || $len > $MAXVALUE ) {
#		print STDERR "ERROR SIZE @ line $.: $_\n!";
#		next;
#	}
	$fraglen[ $len ] ++;
	++ $all;
}
print STDERR "\rDone: $. lines loaded.\n";
close IN;

if( $all == 0 ) {
	print STDERR "ERROR: no valid reads!\n";
	exit 1;
}

#my $maxLen = $#fraglen;
#$maxLen = 600 if $maxLen<600;
my $maxLen = 1000;

## mask the tailing zeros
while( 1 ) {
	last if defined $fraglen[ $maxLen ];
	-- $maxLen;
}
#print STDERR "MAX: $maxLen\n";

print "#length\tcount\tpercent%\tcumulative\n";
my $here = 0;
my $i = 1;
## the commented code will mask the starting zeros
#for( ; $i<=$maxLen; ++$i ) {
#	last if defined $fraglen[$i];
#}
for( ; $i<=$maxLen; ++$i ) {
	$fraglen[$i] = 0 unless defined $fraglen[$i];
	$here += $fraglen[$i];
	print "$i\t", $fraglen[$i], "\t", $fraglen[$i]/$all*100, "\t", $here/$all, "\n";
}

