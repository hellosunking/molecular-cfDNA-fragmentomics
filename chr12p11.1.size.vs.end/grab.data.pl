#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <bed.list> [chr=chr12] [start=34331798] [end=34407798]\n\n";
	exit 2;
}

my $bed   = $ARGV[0];
my $chr   = $ARGV[1] || 'chr12';
my $start = $ARGV[2] || 34331798;
my $end   = $ARGV[3] || 34407798;

my $altchr = $chr;
$altchr =~ s/chr/hs/;	## for handling public datasets

print STDERR "Loading $bed\n";
if( $bed =~ /\.gz$/ ) {
	open BED, "pigz -cd -p 4 $bed | " or die( "$!" );
} else {
	open BED, "$bed" or die( "$!" );
}

while( my $r = <BED> ) {
	chomp($r);
	my @l = split /\t/, $r;
	if( $l[0] eq $chr || $l[0] eq $altchr ) {
		next if $l[2] < $start;
		last if $l[1] > $end;
		print "$chr\t$l[1]\t$l[2]\n";
	}
}

close BED;

