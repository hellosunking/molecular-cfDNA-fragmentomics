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

sid=argv[1];

outfileName = paste0(sid, ".size.vs.end.pdf");
pdf( outfileName, width=10 );

## load data
allU   = read.table(paste(sid,   "all.U.dist", sep="."));
allD   = read.table(paste(sid,   "all.D.dist", sep="."));
longU  = read.table(paste(sid,  "long.U.dist", sep="."));
longD  = read.table(paste(sid,  "long.D.dist", sep="."));
shortU = read.table(paste(sid, "short.U.dist", sep="."));
shortD = read.table(paste(sid, "short.D.dist", sep="."));

## normalize data
allU$V2   = allU$V2   / sum( allU$V2  )*100;
allD$V2   = allD$V2   / sum( allD$V2  )*100;
longU$V2  = longU$V2  / sum( longU$V2 )*100;
longD$V2  = longD$V2  / sum( longD$V2 )*100;
shortU$V2 = shortU$V2 / sum(shortU$V2 )*100;
shortD$V2 = shortD$V2 / sum(shortD$V2 )*100;

## plot 1: all fragments
plot( -1000, xlim=c(-100, 100), ylim=range(allU$V2, allD$V2),
	 xlab="Relative position to nucleosome center (bp)", ylab="Normalized end density",
	 main=paste( sid, ": all fragments", sep=""));
abline( v=c(-95, -73, 73, 95), col='grey', lty=2);
lines(allU$V2~allU$V1, col='grey8',  lwd=2);
lines(allD$V2~allD$V1, col='grey48', lwd=2);
legend( 'top', c('U end', 'D end'), col=c('grey8', 'grey48'),
	   lty=c(1,1), bty='n', cex=1.4, lwd=2 );

## plot 2: long fragments
plot( -1000, xlim=c(-100, 100), ylim=range(longU$V2, longD$V2),
	 xlab="Relative position to nucleosome center (bp)", ylab="Normalized end density",
	 main=paste( sid, ": long fragments (>=170)", sep=""));
abline( v=c(-95, -73, 73, 95), col='grey', lty=2);
lines(longU$V2~longU$V1, col='purple',       lwd=2);
lines(longD$V2~longD$V1, col='springgreen4', lwd=2);
legend( 'top', c('U end', 'D end'), col=c('purple', 'springgreen4'),
	   lty=c(1,1), bty='n', cex=1.4, lwd=2 );

## plot 3: short fragments
plot( -1000, xlim=c(-100, 100), ylim=range(shortU$V2, shortD$V2),
	 xlab="Relative position to nucleosome center (bp)", ylab="Normalized end density",
	 main=paste( sid, ": short fragments (<=147)", sep=""));
abline( v=c(-95, -73, 73, 95), col='grey', lty=2);
lines(shortU$V2~shortU$V1, col='red',  lwd=2);
lines(shortD$V2~shortD$V1, col='blue', lwd=2);
legend( 'top', c('U end', 'D end'), col=c('red', 'blue'),
	   lty=c(1,1), bty='n', cex=1.4, lwd=2 );

## plot 4: merged short vs long
plot( -1000, xlim=c(-100, 100), ylim=range(longU$V2, longD$V2, shortU$V2, shortD$V2),
	 xlab="Relative position to nucleosome center (bp)", ylab="Normalized end density",
	 main=paste( sid, ": short vs long fragments", sep=""));
abline( v=c(-95, -73, 73, 95), col='grey', lty=2);
lines(longU$V2~longU$V1, col='purple',       lwd=2);
lines(longD$V2~longD$V1, col='springgreen4', lwd=2);
lines(shortU$V2~shortU$V1, col='red',  lwd=2);
lines(shortD$V2~shortD$V1, col='blue', lwd=2);
legend( 'top', c('Long - U end', 'Long - D end', 'Short - U end', 'Short - D end'),
	   col=c('purple','springgreen4', 'red', 'blue'), lty=c(1,1,1,1), bty='n', cex=1.4, lwd=2 );

