#
# Author: ahfyth
#
# R script for 
#

argv = commandArgs(T);
if( length(argv) < 1 ) {
	print( 'usage: R --slave --args <sid> <info> < plot.R' );
	q();
}

sid = argv[1];
info= "NA";
if( length(argv) > 1 ) {
	info= argv[2];
}

outfileName=paste(sid, ".size.vs.meth.pdf", sep="");
pdf( outfileName );

high = read.table( paste(sid, ".high.size", sep="") );
low  = read.table( paste(sid, ".low.size", sep="") );

ymax=max( high[,3], low[,3] );

## 0~300 bp, for most cases
plot( low$V3 ~ low$V1, type='l', lwd=3, xlim=c(0, 300), ylim=c(0, ymax), col='blue',
			xlab="Fragment size (bp)", ylab="Frequency (%)",
			main=paste( sid, ", info=", info, sep=""), cex.lab=1.5 );
lines( high$V3 ~ high$V1, lwd=3, col='red' );
legend( 'topright', c('Highly methylated', 'Lowly methylated'), col=c('red', 'blue'), lty=c(1,1), bty='n', cex=1.2 );

#dev.off();
