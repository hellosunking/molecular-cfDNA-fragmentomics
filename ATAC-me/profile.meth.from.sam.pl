#!/usr/bin/perl
#
# Author: Kun Sun @SZBL (sunkun@szbl.ac.cn)
#

use strict;
use warnings;

if( $#ARGV < 1 ) {
	print STDERR "\nUsage: $0 <genome.fa> <in.PE.bam> [in.PE.bam ...]\n\n";
	exit 2;
}

my $genome = shift @ARGV;
my $g = loadGenome( $genome );

my $minQual = 30;

## load SAM files
foreach my $sam ( @ARGV ) {
#	print STDERR "Loading $sam\n";
	if( $sam =~ /bam$/ ) {
		open SAM, "samtools view -@ 4 $sam |" or die( "$!" );
	} else {
		open SAM, "$sam" or die( "$!" );
	}
	while( my $r1 = <SAM> ) {	## PE sam
		my $r2 = <SAM>;
		chomp( $r1 );
		my @read1 = split /\t/, $r1, 12;
		next if $read1[2] eq 'chrM';	## discard chrM reads
		next if $read1[2] =~ /_/;		## discard the chr1_random thing

		next if $read1[4] < $minQual;
		next unless exists $g->{$read1[2]};
		chomp( $r2 );
		my @read2 = split /\t/, $r2, 12;
		#SRR87.sra.14	99	chr20	13909320	42	151M	*	0	0	GATAGTTAATGAAT	JJFJJFAAFJJ	XG:Z:CT
		#SRR87.sra.14	147	chr20	13909354	42	151M	*	0	0	TTTTTTATTATAGT	JJJJJF-FAAA	XG:Z:CT

		## fix seqeunce using CIGAR
		my ($seq1, $qual1) = fix_seq_using_CIGAR( $read1[9], $read1[10], $read1[5] );
		my ($seq2, $qual2) = fix_seq_using_CIGAR( $read2[9], $read2[10], $read2[5] );

		## pool read1 and read2 sequence
		my ($p1, $s1, $q1);
		my ($p2, $s2, $q2);

		if( $read1[3] <= $read2[3] ) {
			($p1, $s1, $q1) = ($read1[3], $seq1, $qual1);
			($p2, $s2, $q2) = ($read2[3], $seq2, $qual2);
		} else {
			($p1, $s1, $q1) = ($read2[3], $seq2, $qual2);
			($p2, $s2, $q2) = ($read1[3], $seq1, $qual1);
		}

		my $merge = '';

		if( $p1 + length($s1) <= $p2 ) {   # no overlap between read1 and read2
			$merge = $s1;
			$merge .= '.' x ($p2-$p1-length($s1));
			$merge .= $s2;
		} else {    # merge read1 and read2
			if( $p1 + length($s1) < $p2+length($s2) ) {
				$merge = substr($s1, 0, $p2-$p1);
				for( my $k=$p2-$p1; $k!=length($s1); ++$k ) {   # overlapped region
					if( substr($q1,$k,1) ge substr($q2,$k+$p1-$p2,1) ) {
						$merge .= substr($s1, $k, 1);
					} else {
						$merge .= substr($s2, $k+$p1-$p2, 1);
					}
				}
				$merge .= substr($s2, length($s1)+$p1-$p2);
			} else {	## R1 completely contains R2, very rare
				$merge = $s1;
			}
		}

		## methylation call
		my $ref = substr( $g->{$read1[2]}, $p1-1, length($merge)+1 );
		my ( $offset, $mSig, $uSig );
		if( $read1[-1] =~ /XG:Z:CT/ ) {	## mapped to watson chain
			( $offset, $mSig, $uSig ) = ( 0, 'C', 'T' );
		} elsif( $read1[-1] =~ /XG:Z:GA/ ) {
			( $offset, $mSig, $uSig ) = ( 1, 'G', 'A' );
		} else {
			print STDERR "Error: Unknown strand ($read1[-1]) for $read1[0]!\n";
			next;
		}
		my $i = -1;
		my ( $CpG, $M, $U, $N ) = (0, 0, 0, 0);
		while( ($i=index($ref, "CG", $i+1)) >= 0 ) {
			++ $CpG;
			my $s = substr( $merge, $i+$offset, 1 );
			if( $s eq $mSig ){ ++ $M; }
			elsif( $s eq $uSig ){ ++ $U; }
			else {
				if( $s eq '.' ) {
					-- $CpG;	## mark this CpG site NO covered
				} else {
					++ $N;
				}
			}
		}
		$p1 --;	## for 0-base BED format
		print join("\t", $read1[2], $p1, $p1+length($merge), $read1[0], $CpG, $M, $U, $N), "\n";
	}
	close SAM;
}

sub fix_seq_using_CIGAR {
	my $seq   = shift;
	my $qual  = shift;
	my $cigar = shift;

	if( $cigar !~ /^\d+M$/ ) {	## there are indels here
		$cigar =~ s/([MID])/$1:/g;
		my @info = split /:/, $cigar;
		my ($mseq, $mqual) = ('', '');	## modified seq and qual
		my $curr = 0;
		foreach my $m ( @info ) {
			if( $m =~ /^(\d+)I$/ ) {	## insertion: CAUTION!!! I WILL DISCARD IT!!!
				$curr += $1;
			} elsif ( $m =~ /^(\d+)D$/ )	{ ## deletion: CAUTION!!! I WILL ADD 'D' to indicate deletion !!!
				$mseq .= 'D' x $1;
				$mqual.= '!' x $1;
			} elsif( $m =~ /^(\d+)M$/ ) {
				$mseq .= substr($seq,  $curr, $1);
				$mqual.= substr($qual, $curr, $1);
				$curr += $1;
			}
		}
		return ($mseq, $mqual);
	} else {	## nothing to update
		return ( $seq, $qual);
	}
}

sub loadGenome {
	my $fasta = shift;

	my %g;
	my $chr = 'NA';
	open IN, "$fasta" or die("$!");
	while( <IN> ) {
		chomp;
		if( /^>(\S+)/ ) {
			$chr = $1;
		} else {
			$g{"$chr"} .= uc $_;
		}
	}
	close IN;

	return \%g;
}

