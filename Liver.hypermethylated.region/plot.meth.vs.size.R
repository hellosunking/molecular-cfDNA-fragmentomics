#
# Author: ahfyth
#
# R script for 
#

argv = commandArgs(T);
if( length(argv) != 1 ) {
	print( 'usage: R --slave --args <sid> < plot.R' );
	q();
}

sid = argv[1];

outfileName=paste0(sid, ".meth.vs.size.pdf");
pdf( outfileName );

hyper = read.table( paste0(sid, ".methed.size") );
hypo  = read.table( paste0(sid, ".unmeth.size") );

## normalization
hyper$V2 = hyper$V2/sum(hyper$V2) * 100;
hypo$V2  = hypo$V2 /sum(hypo$V2)  * 100;

## density
ymax = max( hyper$V2, hypo$V2 );
plot( hyper$V2 ~ hyper$V1, type='l', lwd=2, xlim=c(50, 250), ylim=c(0, ymax), col='red',
			xlab="Fragment size (bp)", ylab="Frequency (%)", main=argv[2] );
lines( hypo$V2 ~ hypo$V1, lwd=2, col='blue' );
legend( 'topleft', c("Hyper", "Hypo"), col=c('red', 'blue'), lty=c(1,1), bty='n', cex=1.1 );

## cumulative density
plot( hyper$V3 ~ hyper$V1, type='l', lwd=2, xlim=c(0, 500), ylim=c(0, 1), col='red',
			xlab="Fragment size (bp)", ylab="Frequency (%)", main=argv[2] );
lines( hypo$V3 ~ hypo$V1, lwd=2, col='blue' );
legend( 'topleft', c("Hyper", "Hypo"), col=c('red', 'blue'), lty=c(1,1), bty='n', cex=1.1 );

