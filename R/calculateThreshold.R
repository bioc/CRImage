calculateThreshold <-
function(imgG,img,method="otsu",window=1,indexWhitePixel=c()){
	message("thresholding222")
	dimensions=dim(imgG)
#exclude white pixels from threshold
	excludeWhite=FALSE 
	if(length(indexWhitePixel)>0){
		message("Exclude white pixel")
		imgG[indexWhitePixel]=-1
		excludeWhite=TRUE
	}
	imgG=array(imgG,dimensions)
	if(method=="otsu"){
		t=localOtsuThreshold(imgG,window,excludeWhite)
	}else if(method=="phansalkar"){
		t=localORThreshold(imgG,img,window,excludeWhite)
	}else{
		message("Thresholding method not known; Otsu thresholding is applied.")
		t=localOtsuThreshold(imgG,img,window,excludeWhite)
	}
	t
}
