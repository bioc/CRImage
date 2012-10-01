segmentImage <-
function(filename="",image=NA,maxShape=NA,minShape=NA,failureRegion=NA,threshold="otsu",numWindows=2,colorCorrection=FALSE,classifyStructures=FALSE,pixelClassifier=NULL,greyscaleImage=0,penClassifier=NULL,referenceHist=NULL){
	if(filename!=""){
		img = readImage(filename)
	}else{
		img=image
	}
	colorMode(img)="color"
	if(length(dim(img))>2){
		imgT=img[,,1]
		whitePixelMask=img[,,1]>0.85 & img[,,2]>0.85 & img[,,3]>0.85
		indexWhitePixel=which(img[,,1]>0.85 & img[,,2]>0.85 & img[,,3]>0.85)
	}else{
		imgT=img
		indexWhitePixel=which(img>0.85)
	}
	numWhitePixel=length(imgT[indexWhitePixel])
	numPixel=length(as.vector(imgT))
	#use these values if defined in order to save time
	imgReduced=NULL
	whitePixelMaskReduced=NULL
	
	#################################################
	if((numWhitePixel/numPixel)<0.9){
		message("Start segmentation.")
		if(classifyStructures==TRUE){
			listStructures=segmentStructures(img,pixelClassifier)
			imgSt=listStructures[[1]]
		}
		if(!is.null(penClassifier)){
			if(is.null(imgReduced)){imgReduced=resize(img,250,250)}
			if(is.null(whitePixelMaskReduced)){whitePixelMaskReduced=imgReduced[,,1]>0.85 & imgReduced[,,2]>0.85 & imgReduced[,,3]>0.85}
			penLabel=classifyPen(imgReduced,penClassifier,whitePixelMaskReduced)
		}
		if(!is.null(referenceHist)){
			if(is.null(imgReduced)){imgReduced=resize(img,250,250)}
			if(is.null(whitePixelMaskReduced)){whitePixelMaskReduced=imgReduced[,,1]>0.85 & imgReduced[,,2]>0.85 & imgReduced[,,3]>0.85}
			distToRef=getImageDistance(imgReduced,referenceHist,whitePixelMaskReduced)
		}
#search large dark regions and color them white
		imgG=channel(img,"grey")
		
		imgB=imgG
		imgB[imgB<0.2]=1
		imgB[imgB !=1]=0
		imgS=imgB
		
		imgS=opening(imgS,makeBrush(5, shape='disc'))
		imgS= bwlabel(imgS)
		imgSd=as.vector(imageData(imgS)+1)
		hc=tabulate(imgSd)
		imgSdN=imageData(imgS)+1
		rm(imgS)
#set the failure region
		if(!is.na(failureRegion)){
			message("Delete failure regions.")
			a=array(hc[imgSdN]<failureRegion,dim(imgSdN))
			imgSdN[a]=-1
		}
		failures=which(imgSdN==-1)
		imgSdN=imgSdN-1
		imgR=img[,,1]
		imgGreen=img[,,2]
		imgBlue=img[,,3]
		failures=which(imgSdN==-1)
		#failure regions
		imgR[failures]=2
		imgGreen[failures]=2
		imgBlue[failures]=2
		img[,,1]=imgR
		img[,,2]=imgGreen
		img[,,3]=imgBlue
		rm(imgR)
		rm(imgBlue)
		rm(imgGreen)
		#failure regions in thre greyscale image are colored white	
		imgG[failures]=1
		imgG[indexWhitePixel]=1
		positivePixels=dim(img[,,1])[1]*dim(img[,,1])[2]-length(indexWhitePixel)-length(failures)
		
		#set the failure region
		#local thresholding, failures are considered as whitePixel
		if(greyscaleImage==0){
			#imgB=calculateThreshold(imgG,img,threshold,window,c(indexWhitePixel,failures))
			#local thresholding, failures are considered as whitePixel
			imgB=createBinaryImage(imgG,img,method=threshold,numWindows=numWindows,whitePixelMask=whitePixelMask)
			writeImage(imgB,"imgB.jpg")
		}else{
			#imgB=calculateThreshold(img[,,greyscaleImage],img,threshold,window,c(indexWhitePixel,failures))
			#local thresholding, failures are considered as whitePixel
			#whitePixelMask[failureMask]=TRUE
			imgB=createBinaryImage(imgG,img,method=threshold,numWindows=numWindows,whitePixelMask=whitePixelMask)
			writeImage(imgB,"imgB.jpg")
		}
#morphological opening to smooth the shape
		imgB=opening(imgB,makeBrush(5, shape='disc'))
		imgS=bwlabel(imgB)
		rm(imgB)
		
		
#segment the image
		imgSdN=imageData(imgS)+1
		numSeq=tabulate(imageData(imgS)+1)
#set minShape
		if(!is.na(minShape)){
			message("Delete min cell nuclei.")
			a=array(numSeq[imgSdN]<minShape,dim(imgSdN))
			imgSdN[a]=1
		}
		imgSdN=imgSdN-1 
#watershed segmentation
		imgT=distmap(imgSdN)
		imgW=watershed(imgT)
		rm(imgT)
		
		if(!is.na(maxShape)){
			message("Resegment large regions.")
			imgSd=as.vector(imageData(imgW)+1)
			hc=tabulate(imgSd)
			rm(imgSd)
			imgSdN=imageData(imgW)+1
			a=array(hc[imgSdN]<maxShape,dim(imgSdN))
			imgSdN[a]=1
			imgSdN=imgSdN-1
			
#segment large segments another time
			imgW[imgSdN >0]=0
			
			imgLS=imgG
			if(greyscaleImage==0){
				imgLS=imgG
			}else{
				imgLS=img[,,greyscaleImage]
			}
			imgLS[imgSdN==0]=0
			imgSdN=as.vector(imgSdN)
			imgLS=as.vector(imgLS)
			imgB_ls=imgSdN
			allSegments=unique(imgSdN)
			allSegments=allSegments[allSegments !=0]
			
			
			sapply(allSegments,function(i){
				   greyPixel=imgLS[imgSdN==i]
				   indexGreyPixel=which(imgSdN==i)
				   localThreshold=calculateOtsu(as.vector(greyPixel[greyPixel !=1]))
				   thresholdedPixel=greyPixel
				   thresholdedPixel[thresholdedPixel<localThreshold]=-1
				   thresholdedPixel[thresholdedPixel!=-1]=0
				   thresholdedPixel[thresholdedPixel==-1]=1
				   imgB_ls[indexGreyPixel]<<-thresholdedPixel
				   })
			imgB_ls=array(imgB_ls,dim(imgG))
			imgB_ls=opening(imgB_ls)
			imgT_ls=distmap(imgB_ls)
			message("Watershed")
			imgW_ls=watershed(imgT_ls)
			rm(imgT_ls)
			maxSegment=max(imgW)
			imgW_ls=imgW_ls+maxSegment
			imgW_ls[imgW_ls==maxSegment]=0
			message("Large segments resegmented.")
			imgW=imgW+imgW_ls
#delete Border cells####################################
			borderLeft=as.vector(unique(imgW[,1]))
			borderRight=as.vector(unique(imgW[,dim(imgW)[2]]))
			borderUp=as.vector(unique(imgW[1,]))
			borderDown=as.vector(unique(imgW[dim(imgW)[1],]))
			borderCells=c(borderLeft,borderRight,borderUp,borderDown)
			cellsToMask=rep(0,max(imgW)+1)
			cellsToMask[borderCells+1]=1
			imgSdBorder=imageData(imgW)+1
			a=array(cellsToMask[imgSdBorder]==1,dim(imgSdBorder))
			imgSdBorder[a]=1
			imgSdBorder=imgSdBorder-1
			imgW=imgSdBorder
########################################################
			
			rm(imgB_ls)
			imgSd=as.vector(imageData(imgW)+1)
			hc=tabulate(imgSd)
			imgSdN=imageData(imgW)+1
			if(!is.na(minShape)){
				a=array(hc[imgSdN]<minShape,dim(imgSdN))
				imgSdN[a]=1
			}
			imgSdN=imgSdN-1
			imgW=imgSdN
			
			
			imgSd=as.vector(imageData(imgW)+1)
			hc=tabulate(imgSd)
			imgSdN=imageData(imgW)+1
			if(!is.na(minShape)){
				a=array(hc[imgSdN]>maxShape,dim(imgSdN))
				imgSdN[a]=1
			}
			imgSdN=imgSdN-1
			imgW=imgSdN
			rm(imgSdN)
		}
		
		message("Relabel regions.")
		regions=unique(as.vector(imgW))
		regions=regions[regions != 0]
		if(length(regions>0)){
#label the nuclei in the right way
			changeLabel=cbind(regions,1:length(regions))
			colnames(changeLabel)=c("oldL","newL")	
			imgWv=as.vector(imgW)
			imgWv_nZ=imgWv[imgWv>0]
			index_imgWv_nZ=which(imgWv>0)
			nV=rep(0,max(imgWv_nZ))
			nV[changeLabel[,1]]=changeLabel[,2]
			imgWv_nZ=nV[imgWv_nZ]
			imgWv[index_imgWv_nZ]=imgWv_nZ
			imgWv=array(imgWv,dim(imgW))
			imgW=imgWv
			rm(imgWv)
#calculate the features
			if(colorCorrection==TRUE){
				meanStdTarget=rbind(c(78.282154,300.371584) ,c(9.694320,-10.856946), c(2.081496,3.614328))
				imgG=colorCorrection(img,meanStdTarget)
			}
			hG=computeFeatures.shape(imgW)
			hF=computeFeatures.shape(imgW)
			allM=computeFeatures.moment(imgW,imgG)
			#zM=computeFeatures.zernikeMoments(imgW,imgG)
			hFgrey=computeFeatures.haralick(imgW,imgG)
			#
			allFeatures=data.frame(stringsAsFactors=FALSE)
			
			allFeatures=cbind(hF,allM,hFgrey)
			allFeatures=allFeatures[allFeatures[,"s.area"]>0,]
			index=1:dim(allFeatures)[1]
			allFeatures=cbind(index,allFeatures)
			cellCoordinates=allFeatures[,c("m.cx","m.cy")]
			cellCoordinates[,1]=as.numeric(format(cellCoordinates[,1],digits=4))
			cellCoordinates[,2]=as.numeric(format(cellCoordinates[,2],digits=4))
			#write.table(allFeatures[,c("g.x","g.y")],"cellCoordinates.txt")
			densityNeighbors=kde2d(cellCoordinates[,"m.cx"],cellCoordinates[,"m.cy"],n=100)
			dV=densityNeighbors$z
			dimnames(dV)=list(densityNeighbors$x,densityNeighbors$y)
			densityValues=interpolate(cellCoordinates,dV)				
			allFeatures=cbind(allFeatures,densityValues)
			
#segment the cytoplasma
			sizeCytoplasma=segmentCytoplasma(img,imgW,indexWhitePixel,imgG,index,allM)
#calculate the number of neighbors of every nuclei
			numberNeighbors=numberOfNeighbors(img,cellCoordinates,allFeatures)
			allFeatures=cbind(allFeatures,numberNeighbors)
			allFeatures=cbind(allFeatures,sizeCytoplasma)
			
			allFeatures=allFeatures[allFeatures[,"s.area"]>0,]
			if(classifyStructures==TRUE){	
#look if nuclei are within structures
				structures=searchStructures(img,imgW,imgSt,index,allM)
				structures=unlist(structures)
				structures=structures[allFeatures[,"s.area"]>0]
				#l=list(img,imgW,allFeatures,indexWhitePixel,TRUE,listStructures,structures)
				l=list()
				l[["image"]]=img
				l[["segmentedImage"]]=imgW
				l[["features"]]=allFeatures
				l[["whitePixelIndex"]]=indexWhitePixel
				l[["classify"]]=TRUE
				l[["listStructures"]]=listStructures
				l[["structures"]]=structures

				if(!is.null(penClassifier)){l[["penLabel"]]=penLabel}
				if(!is.null(referenceHist)){l[["distToRef"]]=distToRef}
			}else{
				l=list()
				l[["image"]]=img
				l[["segmentedImage"]]=imgW
				l[["features"]]=allFeatures
				l[["whitePixelIndex"]]=indexWhitePixel
				l[["classify"]]=TRUE
				if(!is.null(penClassifier)){l[["penLabel"]]=penLabel}
				if(!is.null(referenceHist)){l[["distToRef"]]=distToRef}
				#l=list(img,imgW,allFeatures,indexWhitePixel,TRUE)
			}
			
			message("Segmentation ended.")
			
			l
		}else{
			message("no cells")
			allFeatures=data.frame(stringsAsFactors=FALSE)
			#l=list(img,img,allFeatures,indexWhitePixel,FALSE)
			l=list()
			l[["image"]]=img
			l[["segmentedImage"]]=img
			l[["features"]]=allFeatures
			l[["whitePixelIndex"]]=indexWhitePixel
			l[["classify"]]=FALSE
			if(!is.null(penClassifier)){l[["penLabel"]]=penLabel}
			if(!is.null(referenceHist)){l[["distToRef"]]=distToRef}
			l
		}
	}else{
		message("Almost white. No classicfication applied.")
		allFeatures=data.frame(stringsAsFactors=FALSE)
		#l=list(img,img,allFeatures,indexWhitePixel,FALSE)
		l=list()
		l[["image"]]=img
		l[["segmentedImage"]]=img
		l[["features"]]=allFeatures
		l[["whitePixelIndex"]]=indexWhitePixel
		l[["classify"]]=FALSE
		l[["penLabel"]]="n"
		l
	}
	
}
