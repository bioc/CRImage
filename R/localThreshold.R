localThreshold <-
function(imgG,img){
	imgB=imgG
	numWindows=2
	xWs=round(dim(img)[1]/numWindows)
	yWs=round(dim(img)[2]/numWindows)
	xlw=1
	xrw=round(dim(img)[1]/numWindows)
	yow=1
	yuw=round(dim(img)[2]/numWindows)
	for(iW in 1:numWindows){
		xlw=1
		xrw=round(dim(img)[1]/numWindows)
		for(jW in 1:numWindows){
			if(xrw>dim(img)[1]){xrw=dim(img)[1]}
			imgGL=imgG[xlw:xrw,yow:yuw]
			if(length(as.vector(imgGL[imgGL != 1]))>0){
				globalThreshold=calculateThreshold(as.vector(imgGL[imgGL != 1]))
			}else{
				globalThreshold=0
			}
			imgBL=imgB[xlw:xrw,yow:yuw]
			imgBL[imgGL<globalThreshold]=-1
			imgBL[imgBL != -1]=0
			imgBL[imgBL==-1]=1
			imgB[xlw:xrw,yow:yuw]=imgBL
			rm(imgBL)
			rm(imgGL)
			xlw=xrw+1
			xrw=xrw+xWs
			
		}
		yow=yuw+1
		yuw=yuw+yWs
		if(yuw>dim(img)[2]){yuw=dim(img)[2]}
	}
	imgB
}

