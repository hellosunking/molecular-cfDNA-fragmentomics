#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL
#

use strict;
use warnings;

if( $#ARGV < 1 ) {
	print STDERR "\nUsage: $0 <in.bed> <out.prefix>\n\n";
#	print STDERR "\nThis program is designed to \n\n";
	exit 2;
}

my $SHORT = 147;
my $LONG  = 170;

open Uall, ">$ARGV[1].all.U" or die( "$!" );
open Dall, ">$ARGV[1].all.D" or die( "$!" );

open Ushort, ">$ARGV[1].short.U" or die( "$!" );
open Dshort, ">$ARGV[1].short.D" or die( "$!" );

open Ulong, ">$ARGV[1].long.U" or die( "$!" );
open Dlong, ">$ARGV[1].long.D" or die( "$!" );

#open U152, ">$ARGV[1].152.U" or die( "$!" );
#open D152, ">$ARGV[1].152.D" or die( "$!" );
#open U166, ">$ARGV[1].166.U" or die( "$!" );
#open D166, ">$ARGV[1].166.D" or die( "$!" );

open IN, "$ARGV[0]" or die( "$!" );
while( <IN> ) {
	chomp;
	my @l = split /\t/;	## chr start end extra
	my $size = $l[2] - $l[1];

	$l[1] ++;	## 0-base to 1-base

	print Uall $l[1], "\n";
	print Dall $l[2], "\n";

	if( $size <= $SHORT ) {
		print Ushort $l[1], "\n";
		print Dshort $l[2], "\n";
	} elsif ( $size >= $LONG ) {
		print Ulong $l[1], "\n";
		print Dlong $l[2], "\n";
#	} elsif( $size == 166 ) {
#		print U166 $l[1], "\n";
#		print D166 $l[2], "\n";
#	} elsif( $size == 152 ) {
#		print U152 $l[1], "\n";
#		print D152 $l[2], "\n";
	}
}
close IN;

close Uall;
close Dall;
close Ushort;
close Dshort;
close Ulong;
close Dlong;
#close U166;
#close D166;
#close U152;
#close D152;

