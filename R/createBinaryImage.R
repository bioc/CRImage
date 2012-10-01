createBinaryImage <-
function(imgG,img=NULL,method="otsu",threshold=NULL,numWindows=1,whitePixelMask=c()){
	message("thresholding")
	dimensions=dim(imgG)
#exclude white pixels from threshold
	excludeWhite=FALSE 
	if(length(whitePixelMask)>0){
		message("Exclude white pixel")
		excludeWhite=TRUE
	}
	if(is.null(threshold)){
		imgG=array(imgG,dimensions)
		if(method!="phansalkar"){
			t=localThreshold(imgG=imgG,numWindows=numWindows,excludeWhite=excludeWhite,whitePixelMask=whitePixelMask,method=method)
		}else if(method=="phansalkar"){
			if(is.null(img)){
				stop("Color image not defined.")
			}else{
				t=localORThreshold(imgG=imgG,img=img,numWindows=numWindows,excludeWhite=excludeWhite,whitePixelMask=whitePixelMask)
			}
		}else{
			message("Thresholding method not known; Otsu thresholding is applied.")
			t=localThreshold(imgG=imgG,numWindows=numWindows,excludeWhite=excludeWhite,whitePixelMask=whitePixelMask,method="otsu")
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
