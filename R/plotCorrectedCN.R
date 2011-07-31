plotCorrectedCN <-
function(CN, chr=NULL) {
	if (is.null(chr)) stop("chr must be specified\n")
	par(mfrow=c(2,2), oma=c(0.5, 2, 0.5, 0.5))
	sub.CN <- CN$seg[CN$seg$chrom==chr,]
	x.log <- CN$x[CN$x$Chrom==chr,]
	LRR.range <- range(x.log$Corr.LRR[which(is.finite(x.log$Corr.LRR))])
	plot(x.log[,3] ~ x.log[,2], pch=".", ylim=LRR.range, ylab="LRR", xlab="Position",
		 main="LRR Original Data", col="blue")
	segments(x0=sub.CN$loc.start, x1=sub.CN$loc.end, y0=sub.CN$seg.mean, y1=sub.CN$seg.mean, lwd=3, col="darkgrey")
	plot(x.log[,4] ~ x.log[,2], pch=".", ylab="BAF", xlab="Position",
		 main="BAF Original Data", col="blue")
	plot(x.log[,5] ~ x.log[,2], pch=".", ylim=LRR.range,
		 , ylab="LRR", xlab="Position", main="LRR Corrected Data", col="blue")
	segments(x0=sub.CN$loc.start, x1=sub.CN$loc.end, y0=sub.CN$LRR.tum, y1=sub.CN$LRR.tum, lwd=3, col="darkgrey")
	plot(x.log[,6] ~ x.log[,2], pch=".", ylab="BAF", xlab="Position",
		 main="BAF Corrected Data", col="blue")
	title(paste(CN$seg$ID[1], "Chrom", chr), outer=2)
}

