convertRGBToLAB <-
function(imgT){
	message("convert RGB to LAB")
	imgR=as.vector(imgT[,,1])
	imgG=as.vector(imgT[,,2])
	imgB=as.vector(imgT[,,3])
	channelVectors=rbind(imgR,imgG,imgB)
	rm(imgR)
	rm(imgG)
	rm(imgB)
	channelVectors=channelVectors*255
	rgbToLMS=matrix(c(0.3811,0.1967,0.0241,0.5783,0.7244,0.1288,0.0402,0.0782,0.8444),c(3,3))
	LMS=rgbToLMS %*% channelVectors
	LMS[1,]=LMS[1,]
	LMS[2,]=LMS[2,]
	LMS[3,]=LMS[3,]
	lmsToLab1=matrix(c(1/sqrt(3),0,0,0,1/sqrt(6),0,0,0,1/sqrt(2)),c(3,3))
	lmsToLab2=matrix(c(1,1,1,1,1,-1,1,-2,0),c(3,3))
	
	LAB=lmsToLab1 %*% lmsToLab2 %*% LMS
	rm(LMS)
	imgT[,,1]=array(LAB[1,],dim(imgT[,,1]))
	imgT[,,2]=array(LAB[2,],dim(imgT[,,1]))
	imgT[,,3]=array(LAB[3,],dim(imgT[,,1]))
	rm(LAB)
	imgT
}
