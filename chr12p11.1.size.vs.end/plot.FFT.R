#
# Author: Kun Sun @SZBL (sunkun@szbl.ac.cn)
# Date  :
#
# R script for Fig. 2
#

argv = commandArgs(T);
if( length(argv) != 3 ) {
	print( 'usage: R --slave --args <short.U.dist> <short.D.dist> <out.pdf> < plot.R' );
	q();
}

outfileName=argv[3];
pdf( outfileName );

U = read.table( argv[1] );
U = subset(U, V1 > -100 & V1 < 100 );
lo = loess(U[,2] ~ seq(1,length(U[,2])), span=0.25);
b = U[,2] - predict(lo);
c = spec.pgram(b, pad=0.1, tap=0.1, span=2, plot=F, detrend=T, demean=T)
d = 1/c$freq
cs= c$spec/sum(c$spec);
plot(cs ~ d, type='l', col='red', lwd=2,
		xlim=c(0,50), ylim=c(0, 0.15),
		xlab="Distance (bp)", ylab="Spectrum intensity",
		main="FFT of U/D ends in short-size reads");

D = read.table( argv[2] );
D = subset(D, V1 > -100 & V1 < 100 );
lo = loess(D[,2] ~ seq(1,length(D[,2])), span=0.25);
b = D[,2] - predict(lo);
c = spec.pgram(b, pad=0.1, tap=0.1, span=2, plot=F, detrend=T, demean=T)
d = 1/c$freq
cs= c$spec/sum(c$spec);

lines( cs ~ d, col='blue', lwd=2);
legend( 'topright', c('U', 'D'), col=c('red','blue'), lty=c(1,1), bty='n', cex=1.1 );

