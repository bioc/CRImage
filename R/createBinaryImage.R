createBinaryImage <-
function(imgG,img=NULL,method="otsu",threshold=NA,numWindows=1,whitePixelMask=c()){
	message("thresholding")
	dimensions=dim(imgG)
	#exclude white pixels from threshold
	excludeWhite=FALSE 
	if(length(whitePixelMask)>0){
		message("Exclude white pixel")
		excludeWhite=TRUE
	}
	if(is.na(threshold)){
		imgG=array(imgG,dimensions)
		if(method=="otsu"){
			t=localOtsuThreshold(imgG,numWindows,excludeWhite,whitePixelMask)
		}else if(method=="phansalkar"){
			if(is.null(img)){
				stop("Color image not defined.")
			}else{
				t=localORThreshold(imgG,img,numWindows,excludeWhite)
			}
		}else{
			message("Thresholding method not known; Otsu thresholding is applied.")
			t=localOtsuThreshold(imgG,img,numWindows,excludeWhite)
		}
	}else{
		message("use fixed threshold")
		t=array(0,dim(imgG))
		#find foreground
		foreground=which(imgG<threshold)
		#set foreground to one
		t[foreground]=1
		
	}
	t
	
}

