#!/bin/bash
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :
#
set -o nounset
set -o errexit
#command || { echo "command failed"; exit 1; }

if [ $# -lt 4 ]
then
	echo "Usage: $0 <genome.fasta> <ATAC.peak.bed> <in.nameSrt.sam/bam> <experiment.ID>"
	echo "NOTE: the sam/bam file MUST be sorted by reads name."
	exit 2
fi > /dev/stderr

fasta=$1
peak=$2
bam=$3
sid=$4

perl profile.meth.from.sam.pl $fasta $bam >$sid.meth.bed
bedtools intersect -a $sid.meth.bed -b $peak -v > $sid.non.peak.bed
perl separate.by.meth.pl $sid.non.peak.bed $sid
perl bed2size.pl $sid.high.bed y >$sid.high.size
perl bed2size.pl $sid.low.bed  y >$sid.low.size
R --slave --args $sid < plot.R
