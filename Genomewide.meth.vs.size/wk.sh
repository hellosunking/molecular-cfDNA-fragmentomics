#!/bin/bash
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :
#

set -o nounset
set -o errexit
#command || { echo "command failed"; exit 1; }

if [ $# -lt 3 ]
then
	echo "Usage: $0 <genome.fa> <in.nameSrt.bam> <experiment.ID>"
	echo "NOTE: the input bam file MUST be sorted by read name."
	exit 2
fi > /dev/stderr

fasta=$1
inbam=$2
sid=$3

#g++ -O2 -o call.meth.dir -std=c++11 meth.caller.cpp

mkdir -p $sid
perl pool.and.size.separate.pl $sid $inbam
./call.meth.dir $fasta $sid
perl calc.avg.meth.pl $sid/*.meth.log | sort -k1,1n > $sid.meth
wc -l $sid/*.sam | perl -lane 'if($F[1]=~/(\d+).sam/){print "$1\t$F[0]"}' | sort -k1,1n > $sid.size
R --slave --args $sid < plot.R

