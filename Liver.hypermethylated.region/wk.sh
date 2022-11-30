#!/bin/bash
#
# Author: Kun Sun @SZBL (sunkun@szbl.ac.cn)
#

set -o nounset
set -o errexit
#command || { echo "command failed"; exit 1; }

if [ $# -lt 1 ]
then
	echo "Usage: $0 <in.srt.bam> <experiment.ID>" > /dev/stderr
	exit 2
fi

bam=$1
sid=$2

## 'Liver.hyper.selected.bed' is generated using 'mine.diff.pl',
## which compares Msuite2 output bedgraph files from ENCODE's T-cell and Liver tissue WGBS data

extend=500
while read chr pos1 pos2 extra
do
	let  left=$pos1-$extend
	let right=$pos2+$extend
	samtools view $bam $chr:$left-$right | perl process.sam.pl - $chr:$pos2
done < Liver.hyper.selected.bed >$sid.alleles

perl extract.size.pl $sid.alleles $sid
R --slave --args $sid < plot.meth.vs.size.R

