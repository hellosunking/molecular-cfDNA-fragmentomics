#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
#
# R script for 
#

argv = commandArgs(T);
if( length(argv) != 3 ) {
	print( 'usage: R --slave --args <out.pdf> <in.size> <in.Name> < plot.R' );
	q();
}

outfileName=argv[1];
pdf( outfileName );

dat = read.table( argv[2] );
plot( dat[,2]~dat[,1], type='l', lwd=2, xlim=c(0, 250),
		xlab="Fragment size (bp)", ylab="Frequency (%)", cex.lab=1.2,
		main=paste0("Size distribution for ", argv[3]) );

plot( dat[,2]~dat[,1], type='l', lwd=2, xlim=c(0, 600),
		xlab="Fragment size (bp)", ylab="Frequency (%)", cex.lab=1.2,
		main=paste0("Size distribution for ", argv[3]) );

