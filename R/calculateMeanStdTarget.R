calculateMeanStdTarget <-
function(imgT){
	
	imgR=as.vector(imgT[,,1])
	imgG=as.vector(imgT[,,2])
	imgB=as.vector(imgT[,,3])
	channelVectors=rbind(imgR,imgG,imgB)
	labT=convertRGBToLABOld(channelVectors)
	lSdT=sd(labT[1,])
	lMT=mean(labT[1,])
	
	aSdT=sd(labT[2,])
	aMT=mean(labT[2,])
	
	bSdT=sd(labT[3,])
	bMT=mean(labT[3,])
	
	rbind(c(lSdT,lMT),c(aSdT,aMT),c(bSdT,bMT))
	
}
