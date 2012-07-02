localOtsuThresholdNew <-
function(imgG,numWindows,excludeWhite=FALSE, whitePixelMask){

	imgB=imgG
	numWindows=numWindows
	xWs=round(dim(imgG)[1]/numWindows)
	yWs=round(dim(imgG)[2]/numWindows)
	xlw=1
	xrw=round(dim(imgG)[1]/numWindows)
	yow=1
	yuw=round(dim(imgG)[2]/numWindows)
	for(iW in 1:numWindows){
		xlw=1
		xrw=round(dim(imgG)[1]/numWindows)
		for(jW in 1:numWindows){
			if(xrw>dim(imgG)[1]){xrw=dim(imgG)[1]}
			imgGL=imgG[xlw:xrw,yow:yuw]
			# shall white pixel be excluded from thresholding
			if(excludeWhite){
				whiteMaskGL=whitePixelMask[xlw:xrw,yow:yuw]
			}
			if(excludeWhite){
				if(length(as.vector(imgGL[whiteMaskGL == TRUE]))>0){
							globalThreshold=calculateOtsu(as.vector(imgGL[whiteMaskGL == TRUE]))
				}else{
							globalThreshold=NA
				}
			}else{
				if(length(as.vector(imgGL))>0){
							globalThreshold=calculateOtsu(as.vector(imgGL))
				}else{
							globalThreshold=NA
				}
			}
			
			if(is.nan(globalThreshold)){
				imgB[xlw:xrw,yow:yuw]=0
			}else{
				imgBL=array(0,dim(imgGL))
				#find foreground
				foreground=which(imgGL<globalThreshold)
				#set foreground to one
				imgBL[foreground]=1
				#set white pixel or failures to background
				imgBL[imgGL== -1]=0
				
				imgB[xlw:xrw,yow:yuw]=imgBL
				rm(imgBL)
				rm(imgGL)
			}
				xlw=xrw+1
				xrw=xrw+xWs
			
			
		}
		yow=yuw+1
		yuw=yuw+yWs
		if(yuw>dim(imgG)[2]){yuw=dim(imgG)[2]}
	}
	imgB
}

