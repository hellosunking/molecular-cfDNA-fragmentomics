#!/bin/bash
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# This program is for paragraph 1 in Results.
#
# Requirements: R, bedtools
#

set -o nounset
set -o errexit
#command || { echo "command failed"; exit 1; }

if [ $# -lt 2 ]
then
	echo "Usage: $0 <experiment.ID> <bed.list>"
	echo "experiment.ID is the ID of your experiment(s), and the main output files will have this prefix."
	echo "bed.list is a file containing the sampleID and path to bed file. 'bed.list.example' is an example."
	exit 2
fi > /dev/stderr

sid=$1
bedlist=$2

## nucleosome track: DANPOS using MNASE-seq data on Lymphoblastoid cell line GM18522 (Gaffney PLOS Genetics 2012)
if [ ! -s Genomewide.nucleosomes.hg38.bed ]
then
	wget https://download.cncb.ac.cn/nucmap/organisms/v1/Homo_sapiens/byDataType/Nucleosome_peaks_DANPOS/Homo_sapiens.hsNuc0390101.nucleosome.DANPOSPeak.bed.gz
	zcat Homo_sapiens.hsNuc0390101.nucleosome.DANPOSPeak.bed.gz | perl -lane 'next unless $F[0]=~/chr\d+$/; $c=($F[1]+$F[2])>>1; print join("\t", $F[0], $c-101, $c+100, "DANPOS")' | sort -k1,1 -k2,2n >Genomewide.nucleosomes.hg38.bed
fi

## long/short/exact size reads
echo "Overlap Reads vs Nucleosomes ..."
while read sname bedfile
do
	bedtools intersect -a $bedfile -b Genomewide.nucleosomes.hg38.bed -sorted -wao | perl anno.end.pl $sname - &
done < $2
wait

echo "Collecting data ..."
perl pool.end.pl $bedlist $sid

## plot the ends
R --slave --args $sid < plot.end.R

