#
# Author: ahfyth
#
# R script for 
#

argv = commandArgs(T);
if( length(argv) != 1 ) {
	print( 'usage: R --slave --args <prefix> < plot.R' );
	q();
}

outfileName = paste0( argv[1], ".chr12p11.1.meth.vs.end.pdf" );
pdf( outfileName, width=10 );

hyperU = read.table(paste0(argv[1], "hyper.Uend.dist"));
hyperD = read.table(paste0(argv[1], "hyper.Dend.dist"));
hypoU  = read.table(paste0(argv[1], "hypo.Uend.dist" ));
hypoD  = read.table(paste0(argv[1], "hypo.Dend.dist" ));

hyperU$V2 = hyperU$V2 / sum(hyperU$V2)*100;
hyperD$V2 = hyperD$V2 / sum(hyperD$V2)*100;
hypoU$V2  = hypoU$V2  / sum(hypoU$V2 )*100;
hypoD$V2  = hypoD$V2  / sum(hypoD$V2 )*100;

## plot 1: hypermethylated
plot( -1000, xlim=c(-100, 100), ylim=c(0, max(hyperU$V2, hyperD$V2)),
	 xlab="Relative position to nucleosome center (bp)", ylab="Normalized end density",
	 main="Fragments with hypermethylated CpG sites");
abline( v=c(-95, -73, 73, 95), col='grey', lty=2);
lines(hyperU$V2~hyperU$V1, col='purple',       lwd=2);
lines(hyperD$V2~hyperD$V1, col='springgreen4', lwd=2);
legend( 'top', c('Hyper - U end', 'Hyper - D end'), col=c('purple','springgreen4'),
	   lty=c(1,1), bty='n', cex=1.2, lwd=2 );

## plot 2: hypomethylated
plot( -1000, xlim=c(-100, 100), ylim=c(0, max(hypoU$V2, hypoD$V2)),
	 xlab="Relative position to nucleosome center (bp)", ylab="Normalized end density",
	 main="Fragments with hypomethylated CpG sites");
abline( v=c(-95, -73, 73, 95), col='grey', lty=2);
lines(hypoU$V2~hypoU$V1, col='red',  lwd=2);
lines(hypoD$V2~hypoD$V1, col='blue', lwd=2);
legend( 'top', c('Hypo - U end', 'Hypo - D end'), col=c('red', 'blue'),
	   lty=c(1,1), bty='n', cex=1.2, lwd=2 );

