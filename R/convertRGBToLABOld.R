convertRGBToLABOld <-
function(channelVectors){
	channelVectors=channelVectors*255
	rgbToLMS=matrix(c(0.3811,0.1967,0.0241,0.5783,0.7244,0.1288,0.0402,0.0782,0.8444),c(3,3))
	LMS=rgbToLMS %*% channelVectors
	LMS[1,]=LMS[1,]
	LMS[2,]=LMS[2,]
	LMS[3,]=LMS[3,]
	lmsToLab1=matrix(c(1/sqrt(3),0,0,0,1/sqrt(6),0,0,0,1/sqrt(2)),c(3,3))
	lmsToLab2=matrix(c(1,1,1,1,1,-1,1,-2,0),c(3,3))
	
	LAB=lmsToLab1 %*% lmsToLab2 %*% LMS
	LAB
}
