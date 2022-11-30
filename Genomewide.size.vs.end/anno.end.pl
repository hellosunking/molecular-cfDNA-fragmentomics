#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL
#

use strict;
use warnings;

if( $#ARGV < 1 ) {
	print STDERR "\nUsage: $0 <sid> <in.ol> [in.ol...]\n\n";
#	print STDERR "\nThis program is designed to \n\n";
	exit 2;
}

my $prefix = shift;

my $IMPOSSIBLE = 1e9;
my $SHORT = 147;
my $LONG  = 170;
my $minQ  = 20;	## min MAPQ

my (%Ua, %Us, %Ul);
my (%Da, %Ds, %Dl);

foreach my $file ( @ARGV ) {
	open IN, "$file" or die( "$!" );
	while( <IN> ) {
		my @l = split /\t/;	##chr1	10016	10109	30	-	chr1	10005	10206	DANPOS	93
		next if $l[3] < $minQ || $l[5] eq '.';

		my $nc = $l[7] - 100;	## nucleosome center
		my ($U, $D) = ($IMPOSSIBLE, $IMPOSSIBLE);

		if( $l[1] >= $l[6] ) {	## left hand is within the nucleosome
			$U = $l[1] + 1 - $nc;
		}
		if( $l[2] <= $l[7] ) {	## right hand is within the nucleosome
			$D = $l[2] - $nc;
		}

		my $size = $l[2] - $l[1];

		if( $U != $IMPOSSIBLE ) {
			$Ua{$U} ++;
			$Us{$U} ++ if $size <= $SHORT;
			$Ul{$U} ++ if $size >= $LONG;
		}

		if( $D != $IMPOSSIBLE ) {
			$Da{$D} ++;
			$Ds{$D} ++ if $size <= $SHORT;
			$Dl{$D} ++ if $size >= $LONG;
		}

	}
	close IN;
}

open Uall,   ">$prefix.all.U.dist"   or die( "$!" );
open Dall,   ">$prefix.all.D.dist"   or die( "$!" );
open Ushort, ">$prefix.short.U.dist" or die( "$!" );
open Dshort, ">$prefix.short.D.dist" or die( "$!" );
open Ulong,  ">$prefix.long.U.dist"  or die( "$!" );
open Dlong,  ">$prefix.long.D.dist"  or die( "$!" );

foreach my $i (-100..100) {
	print Uall   join("\t", $i, $Ua{$i}||0), "\n";
	print Dall   join("\t", $i, $Da{$i}||0), "\n";
	print Ushort join("\t", $i, $Us{$i}||0), "\n";
	print Dshort join("\t", $i, $Ds{$i}||0), "\n";
	print Ulong  join("\t", $i, $Ul{$i}||0), "\n";
	print Dlong  join("\t", $i, $Dl{$i}||0), "\n";
}

close Uall;
close Dall;
close Ushort;
close Dshort;
close Ulong;
close Dlong;


