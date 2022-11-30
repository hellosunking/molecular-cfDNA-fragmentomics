#
# Author: Kun Sun @SZBL (sunkun@szbl.ac.cn)
# Date  :
#
# R script for 
#

argv = commandArgs(T);
if( length(argv) != 1 ) {
	print( 'usage: R --slave --args <sid> < plot.R' );
	q();
}

sid=argv[1];
outfileName=paste0(sid, '.meth.vs.size.pdf');

pdf( outfileName, width=8, height=6 );
par( mar = c(5, 4, 3, 4)+0.3 );

meth = read.table( paste0( sid, ".meth") );
size = read.table( paste0( sid, ".size") );
size$V2 = size$V2 / sum( size$V2 ) * 100;

plot( meth$V4 ~ meth$V1, type='l', lwd=3, col='red', xlim=c(50,200), #axes=F,
		xlab="Fragment size (bp)", ylab="", cex.lab=1.5, main=sid );
axis(side=2, at=pretty(range(meth$V4)), col='red', col.ticks='red');
mtext("Methylation level (%)", side=2, line=3, col="red", cex=1.5);

par( new=T );
plot( size$V2 ~ size$V1, type='l', lwd=3, col='blue', xlim=c(50,200), axes=F,
		xlab="", ylab="" );
axis(side=4, at=pretty(range(size$V2)), col='blue', col.ticks='blue');
mtext("Size frequency (%)", side=4, line=3,  col="blue", cex=1.5);

