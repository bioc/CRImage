colorCorrectionOld <-
function(imgO,meanStdTarget){
	message("Color correction old.")
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
	indexWhitePixel=which(imgO[,,1]>0.85 & imgO[,,2]>0.85 & imgO[,,3]>0.85)
#do not transform white pixel
	imgR=imgR[-indexWhitePixel]
	imgG=imgG[-indexWhitePixel]
	imgB=imgB[-indexWhitePixel]
	
	channelVectors=rbind(imgR,imgG,imgB)
	
#find nonwhite pixel
	nw_index=which(!(imgO[,,1]>0.85 & imgO[,,2]>0.85 & imgO[,,3]>0.85))
#convert in the LAB color space
	labO=convertRGBToLABOld(channelVectors)
#calculate mean and standard deviation for the differen channels
	lSdO=sd(labO[1,])
	lMO=mean(labO[1,])
	
	aSdO=sd(labO[2,])
	aMO=mean(labO[2,])
	
	bSdO=sd(labO[3,])
	bMO=mean(labO[3,])
	
	
	
#standardize
	labOm=labO
	labOm[1,]=((labOm[1,]-lMO)/lSdO)*lSdT + lMT
	labOm[2,]=((labOm[2,]-aMO)/aSdO)*aSdT + aMT
	labOm[3,]=((labOm[3,]-bMO)/bSdO)*bSdT + bMT
	
#convert back
	RGB=convertLABToRGBOld(labOm)
	RGB=RGB/255
	imgON=imgO
	imgR=as.vector(imgO[,,1])
	imgG=as.vector(imgO[,,2])
	imgB=as.vector(imgO[,,3])
	if(length(nw_index)>0){
		if(length(RGB[1,])>0){
			imgR[nw_index]=RGB[1,]
		}
		if(length(RGB[2,])>0){
			imgG[nw_index]=RGB[2,]
		}
		if(length(RGB[3,])>0){
			imgB[nw_index]=RGB[3,]
		}
	}
	imgON[,,1]=array(imgR,dim(imgO)[1:2])
	imgON[,,2]=array(imgG,dim(imgO)[1:2])
	imgON[,,3]=array(imgB,dim(imgO)[1:2])
	l=list(imgON,indexWhitePixel)
}
