#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL
#

use strict;
use warnings;

if( $#ARGV < 1 ) {
	print STDERR "\nUsage: $0 <out.prefix> <in.bed.gz> [in.bed.gz]\n\n";
#	print STDERR "\nThis program is designed to \n\n";
	exit 2;
}

my $SHORT = 147;
my $LONG  = 170;

my $prefix = shift;

open Uall, ">>$prefix.all.U.bed" or die( "$!" );
open Dall, ">>$prefix.all.D.bed" or die( "$!" );

open Ushort, ">>$prefix.short.U.bed" or die( "$!" );
open Dshort, ">>$prefix.short.D.bed" or die( "$!" );

open Ulong, ">>$prefix.long.U.bed" or die( "$!" );
open Dlong, ">>$prefix.long.D.bed" or die( "$!" );

#open U152, ">$ARGV[1].152.U" or die( "$!" );
#open D152, ">$ARGV[1].152.D" or die( "$!" );
#open U166, ">$ARGV[1].166.U" or die( "$!" );
#open D166, ">$ARGV[1].166.D" or die( "$!" );

foreach my $file ( @ARGV ) {
	open IN, "less $file |" or die( "$!" );
	while( <IN> ) {
		chomp;
		my @l = split /\t/;	## chr start end extra
		my $size = $l[2] - $l[1];

		## mind the 0-base and 1-base thing
		print Uall join($l[0], $l[1], $l[1]+1), "\n";
		print Dall join($l[0], $l[2]-1, $l[2]), "\n";

		if( $size <= $SHORT ) {
			print Ushort join($l[0], $l[1], $l[1]+1), "\n";
			print Dshort join($l[0], $l[2]-1, $l[2]), "\n";
		} elsif ( $size >= $LONG ) {
			print Ulong join($l[0], $l[1], $l[1]+1), "\n";
			print Dlong join($l[0], $l[2]-1, $l[2]), "\n";
		}
	}
	close IN;
}

close Uall;
close Dall;
close Ushort;
close Dshort;
close Ulong;
close Dlong;

