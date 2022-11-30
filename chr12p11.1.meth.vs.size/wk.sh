#!/bin/bash
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
#
#

set -o nounset
set -o errexit
#command || { echo "command failed"; exit 1; }

if [ $# -lt 4 ]
then
	echo "Usage: $0 <in.Msuite.bam> <min.CpG.number> <out.prefix> <chr12.fasta>" > /dev/stderr
	exit 2
fi

inbam=$1
minCpG=$2
prefix=$3
chr12fasta=$4

#echo Use CpG No. cutoff: $minCpG
## the nucleosome centers in chr12p11.1 is obtained using the following command:
#zcat Homo_sapiens.hsNuc0390101.nucleosome.DANPOSPeak.bed.gz | perl -lane 'next unless $F[0] eq "chr12" && $F[1]>=34331798 && $F[2]<=34407798; $c=($F[1]+$F[2])>>1; print join("\t", "chr12", $c-1, $c, "DANPOS")' >chr12p11.1.nucleosome.center.bed
nuctrack=chr12p11.1.nucleosome.center.bed
chr=chr12
sPos=34331000
ePos=34408000

## grab data
samtools view $inbam $chr:$sPos-$ePos | perl filter.PE.sam.pl - > $prefix.within.chr12p11.1.sam

## extract fragment end and methylation information
perl Msuite.sam2bed.PE.pl $prefix.within.chr12p11.1.sam $prefix.result $chr12fasta $minCpG

#### split the data based on methylation level
#### Hyper>=0.8, Hypo<=0.2
cat $prefix.result | perl -lane 'print $F[1] if $F[4]/$F[3]>=0.8' >$prefix.hyper.Uend
cat $prefix.result | perl -lane 'print $F[2] if $F[4]/$F[3]>=0.8' >$prefix.hyper.Dend
cat $prefix.result | perl -lane 'print $F[1] if $F[4]/$F[3]<=0.2' >$prefix.hypo.Uend
cat $prefix.result | perl -lane 'print $F[2] if $F[4]/$F[3]<=0.2' >$prefix.hypo.Dend

## get BED and size
cat $prefix.result | perl -lane 'print "chr12\t$F[1]\t$F[2]" if $F[4]/$F[3]>=0.8' | perl bed2size.pl - >$prefix.hyper.size
cat $prefix.result | perl -lane 'print "chr12\t$F[1]\t$F[2]" if $F[4]/$F[3]<=0.2' | perl bed2size.pl - >$prefix.hypo.size

## annotate the fragment ends using nucleosome track
perl anno.end.pl $nuctrack $prefix.hyper.Uend >$prefix.hyper.Uend.anno 2>$prefix.hyper.Uend.dist
perl anno.end.pl $nuctrack $prefix.hyper.Dend >$prefix.hyper.Dend.anno 2>$prefix.hyper.Dend.dist
perl anno.end.pl $nuctrack $prefix.hypo.Uend  >$prefix.hypo.Uend.anno  2>$prefix.hypo.Uend.dist
perl anno.end.pl $nuctrack $prefix.hypo.Dend  >$prefix.hypo.Dend.anno  2>$prefix.hypo.Dend.dist

## plot
R --slave --args $prefix < plot.ends.R
R --slave --args $prefix.chr12p11.1.meth.vs.size.pdf $prefix $prefix.hyper.size $prefix.hypo.size < plot.size.R

