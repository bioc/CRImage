colorCorrection <-
function(imgO,meanStdTarget,whiteMask=c()){
	message("Color correction.")
#standard deviation of the l channel 
	lSdT=meanStdTarget[1,1]
# mean of the l-space
	lMT=meanStdTarget[1,2]
	
#mean and standard deviation of the a channel
	aSdT=meanStdTarget[2,1]
	aMT=meanStdTarget[2,2]
#mean and standard deviation of the b channel
	bSdT=meanStdTarget[3,1]
	bMT=meanStdTarget[3,2]
	
#seperate the different channels
	imgR=as.vector(imgO[,,1])
	imgG=as.vector(imgO[,,2])
	imgB=as.vector(imgO[,,3])
#find white pixels
	if(length(whiteMask)>0){
		nwMask=!whiteMask
	}else{
		nwMask=matrix(TRUE,dim(imgO[,,1])[1],dim(imgO[,,1])[2])
	}
	indexWhitePixel=which(nwMask==FALSE)
	
	channelVectors=rbind(imgR,imgG,imgB)
	rm(imgR)
	rm(imgB)
	rm(imgG)
	
#convert in the LAB color space
	labO=convertRGBToLAB(imgO)
#calculate mean and standard deviation for the differen channels
	labO_l=as.vector(labO[,,1])
	labO_a=as.vector(labO[,,2])
	labO_b=as.vector(labO[,,3])
	if(length(indexWhitePixel)>0){
		labO_l=labO_l[-indexWhitePixel]
		labO_a=labO_a[-indexWhitePixel]
		labO_b=labO_b[-indexWhitePixel]
	}
	lSdO=sd(labO_l)
	lMO=mean(labO_l)
	
	aSdO=sd(labO_a)
	aMO=mean(labO_a)
	
	bSdO=sd(labO_b)
	bMO=mean(labO_b)
	
	rm(labO_l)
	rm(labO_a)
	rm(labO_b)
	
#standardize
	labO[,,1]=((labO[,,1]-lMO)/lSdO)*lSdT + lMT
	labO[,,2]=((labO[,,2]-aMO)/aSdO)*aSdT + aMT
	labO[,,3]=((labO[,,3]-bMO)/bSdO)*bSdT + bMT
	
#convert back
	RGB=convertLABToRGB(labO)
	rm(labO)
	RGB_r=RGB[,,1]
	RGB_g=RGB[,,2]
	RGB_b=RGB[,,3]
	
	imgR=imgO[,,1]
	imgG=imgO[,,2]
	imgB=imgO[,,3]
	imgR[nwMask]=RGB_r[nwMask]
	imgG[nwMask]=RGB_g[nwMask]
	imgB[nwMask]=RGB_b[nwMask]
	imgO[,,1]=imgR
	imgO[,,2]=imgG
	imgO[,,3]=imgB
	imgO
}
