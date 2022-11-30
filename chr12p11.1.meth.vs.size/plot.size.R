#
# Author: ahfyth
#
# R script for 
#

argv = commandArgs(T);
if( length(argv) != 4 ) {
	print( 'usage: R --slave --args <out.pdf> <main> <hyper.size> <hypo.size> < plot.R' );
	q();
}

outfileName=argv[1];
pdf( outfileName );

hyper = read.table( argv[3] );
hypo  = read.table( argv[4] );

ymax = max( hyper[,3], hypo[,3] );

## density
plot( hyper[,3] ~ hyper[,1], type='l', lwd=2, xlim=c(0, 250), ylim=c(0, ymax), col='red',
			xlab="Fragment size (bp)", ylab="Frequency (%)", main=argv[2] );
lines( hypo[,3]~hypo[,1], lwd=2, col='blue' );
legend( 'topleft', c("Hyper", "Hypo"), col=c('red', 'blue'), lty=c(1,1), bty='n', cex=1.1 );

## cumulative freq
plot( hyper[,4] ~ hyper[,1], type='l', lwd=2, xlim=c(0, 250), ylim=c(0, 1), col='red',
			xlab="Fragment size (bp)", ylab="Frequency (%)", main=argv[2] );
lines( hypo[,4]~hypo[,1], lwd=2, col='blue' );
legend( 'topleft', c("Hyper", "Hypo"), col=c('red', 'blue'), lty=c(1,1), bty='n', cex=1.1 );

