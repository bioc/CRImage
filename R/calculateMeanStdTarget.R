calculateMeanStdTarget <-
function(imgT){
	labT=convertRGBToLAB(imgT)
	lSdT=sd(as.vector(labT[,,1]))
	lMT=mean(as.vector(labT[,,1]))
	aSdT=sd(as.vector(labT[,,2]))
	aMT=mean(as.vector(labT[,,2]))
	bSdT=sd(as.vector(labT[,,3]))
	bMT=mean(as.vector(labT[,,3]))
	rbind(c(lSdT,lMT),c(aSdT,aMT),c(bSdT,bMT))
	
}
