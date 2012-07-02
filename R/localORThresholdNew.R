localORThresholdNew <-
function(imgG,img,numWindows,excludeWhite=FALSE,whitePixelMask){
	message("Phansalkar threshold")
	imgR=as.vector(img[,,1])
	imgGreen=as.vector(img[,,2])
	imgB=as.vector(img[,,3])
	channelVectors=rbind(imgR,imgGreen,imgB)
	labV=convertRGBToLAB(img)
	imgLAB=sqrt(as.vector(labV[,,1])^2+as.vector(labV[,,2])^2+as.vector(labV[,,3])^2)
	imgLAB=matrix(imgLAB/max(imgLAB),dim(imgG))
	imgB=imgG
	imgB_rgb=imgG
	imgB_lab=imgLAB
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
			imgLABL=imgLAB[xlw:xrw,yow:yuw]
			if(excludeWhite){
				whiteMaskGL=whitePixelMask[xlw:xrw,yow:yuw]
			}
			if(excludeWhite){
				#exclude white pixel from thresholding
				if(length(as.vector(imgGL[whiteMaskGL == TRUE]))>0){
					iccsp_rgb=Phansalkar_threshold(as.vector(imgGL[whiteMaskGL == TRUE]))
					iccsp_lab=Phansalkar_threshold(as.vector(imgLABL[whiteMaskGL == TRUE]))
				}else{
					iccsp_rgb=NA
					iccsp_lab=NA
				}
			}else{
				if(length(as.vector(imgGL))>0){
					iccsp_rgb=Phansalkar_threshold(as.vector(imgGL))
					iccsp_lab=Phansalkar_threshold(as.vector(imgLABL))
				}else{
					iccsp_rgb=NA
					iccsp_lab=NA
				}
			}
			
			if(is.nan(iccsp_lab)&&is.nan(iccsp_rgb)){
				imgB[xlw:xrw,yow:yuw]=0
			}else{
				imgBL_rgb=imgB_rgb[xlw:xrw,yow:yuw]
				imgBL_rgb[imgGL<iccsp_rgb]=-1
				imgBL_rgb[imgBL_rgb != -1]=0
				imgBL_rgb[imgBL_rgb==-1]=1
			
				imgBL_lab=imgB_lab[xlw:xrw,yow:yuw]
				imgBL_lab[imgLABL<iccsp_lab]=-1
				imgBL_lab[imgBL_lab != -1]=0
				imgBL_lab[imgBL_lab==-1]=1
			
				sizeImgBL_rgb=dim(imgBL_rgb)
				imgB_temp=array(0,c(sizeImgBL_rgb[2],sizeImgBL_rgb[1]))
				imgB_temp[imageData(imgBL_rgb)==1]=1
				imgB_temp[imageData(imgBL_lab)==1]=1
				imgB[xlw:xrw,yow:yuw]=imgB_temp
				rm(imgB_temp)
				rm(imgGL)
			}
				xlw=xrw+1
				xrw=xrw+xWs
			
			
		}
		yow=yuw+1
		yuw=yuw+yWs
		if(yuw>dim(img)[2]){yuw=dim(img)[2]}
	}
	imgB
}

