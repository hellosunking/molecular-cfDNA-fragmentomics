#!/bin/bash
#
# Author: Kun Sun @SZBL (sunkun@szbl.ac.cn)
#

set -o nounset
set -o errexit
#command || { echo "command failed"; exit 1; }

if [ $# -lt 2 ]
then
	echo "Usage: $0 <in.bed> <experiment.ID>" > /dev/stderr
	exit 2
fi

bed=$1
sid=$2

#### nucleosome center loci
chr=chr12
sPos=34331000
ePos=34408000

## DANPOS using MNASE-seq data on Lymphoblastoid cell line GM18522 (Gaffney PLOS Genetics 2012)
## the nucleosome centers in chr12p11.1 is obtained using the following command:
#zcat Homo_sapiens.hsNuc0390101.nucleosome.DANPOSPeak.bed.gz | perl -lane 'next unless $F[0] eq "chr12" && $F[1]>=34331798 && $F[2]<=34407798; $c=($F[1]+$F[2])>>1; print join("\t", "chr12", $c-1, $c, "DANPOS")' >chr12p11.1.nucleosome.center.bed

target=chr12p11.1.nucleosome.center.bed

## grab reads in chr12p11.1
perl grab.data.pl $bed $chr $sPos $ePos > $sid.within.chr12p11.1.bed

## long/short/exact size reads
perl extract.end.based.on.size.pl $sid.within.chr12p11.1.bed $sid

## annotate the fragment ends using nucleosome track
for C in short long
do
	perl anno.end.pl $target $sid.$C.U >$sid.$C.U.anno 2>$sid.$C.U.dist
	perl anno.end.pl $target $sid.$C.D >$sid.$C.D.anno 2>$sid.$C.D.dist
done

## plot the ends
R --slave --args $sid < plot.end.R

## 'plot.FFT.R' is called manually on pooled control samples

