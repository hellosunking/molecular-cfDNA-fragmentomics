#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;
#use KSLIB::loadGenome qw/loadGenome/;
#use KSLIB::cigarUtil qw/fix_seq_from_CIGAR extract_size_from_CIGAR/;
#use KSLIB::Digitalize qw/digitalize/;

if( $#ARGV < 1) {
	print STDERR "\nUsage: $0 <bed.list> <sid>\n\n";
#	print STDERR "\nThis program is designed to \n\n";
	exit 2;
}

our (%Ua, %Da, %Us, %Ds, %Ul, %Dl);

open IN, "$ARGV[0]" or die( "$!" );
while( <IN> ) {
	chomp;
	my @l = split /\t/;	## sample bed.path
	load_dist( "$l[0].all.U.dist",   \%Ua );
	load_dist( "$l[0].all.D.dist",   \%Da );
	load_dist( "$l[0].short.U.dist", \%Us );
	load_dist( "$l[0].short.D.dist", \%Ds );
	load_dist( "$l[0].long.U.dist",  \%Ul );
	load_dist( "$l[0].long.D.dist",  \%Dl );
}
close IN;

write_dist( "$ARGV[1].all.U.dist",   \%Ua );
write_dist( "$ARGV[1].all.D.dist",   \%Da );
write_dist( "$ARGV[1].short.U.dist", \%Us );
write_dist( "$ARGV[1].short.D.dist", \%Ds );
write_dist( "$ARGV[1].long.U.dist",  \%Ul );
write_dist( "$ARGV[1].long.D.dist",  \%Dl );


sub load_dist {
	my $file = shift;
	my $end  = shift;

	open DIST, "$file" or die( "$!" );
	while( <DIST> ) {
		chomp;
		my @l = split /\t/;
		$end->{$l[0]} += $l[1];
	}
	close DIST;
}

sub write_dist {
	my $file = shift;
	my $end  = shift;

	open DIST, ">$file" or die( "$!" );
	foreach my $i (-100..100) {
		print DIST join("\t", $i, $end->{$i}||0), "\n";
	}
	close DIST;
}

