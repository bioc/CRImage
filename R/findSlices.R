findSlices <-
function(imgFolder,pathToOutputFolder,numSlides,fontSize=10){
	smallImage=readImage(file.path(imgFolder,"SlideThumb.jpg"))
#if not in color mode, convert to color mode
	
#segment the thumbnail to find the sections in the image
	pathToImgFolder=imgFolder
	imgG=channel(smallImage,"gray")
	globalThreshold=calculateOtsu(as.vector(imgG))
	imgG[imgG<globalThreshold]=-1
	imgG[imgG !=-1]=0
	imgG[imgG==-1]=1
	imgG=fillHull(imgG)
	imgG=closing(imgG,makeBrush(10,shape="diamond"))
	imgG=opening(imgG,makeBrush(10,shape="diamond"))
	imgS=bwlabel(imgG)
	hF=hullFeatures(imgS)[,c("g.x","g.y")]
#if the image is found two times (should not happen)
	if(length(numSlides)>1){
		numSlides=numSlides[1]
	}
	
	if(numSlides==1){
		if(is.null(dim((hF)))){
			centers=matrix(c(1,1))
		}else{
			centers=matrix(c(mean(hF[,1])),c(1))
		}
	}else{
#cluster the segments to the number of sections
		if(nrow(hF)>numSlides){
			"centers=kmeans(hF[,1],as.numeric(numSlides),iter.max=2000)$centers"
			centersIndex=cutree(hclust(dist(hF[,1])),k=as.numeric(numSlides))
			centers=matrix(0,as.numeric(numSlides))
			for(i in 1:as.numeric(numSlides)){
				centers[i]=mean(as.numeric(as.character(hF[centersIndex==i,1])))
			}
		}else{
			if(nrow(hF) == numSlides){
				centers=matrix(as.numeric(as.character(hF[,1])),as.numeric(numSlides))
			}
		}
	}
	sliceColors=col2rgb(c("white","green","blue","red","green","yellow","orange","black","brown","purple4","darkolivegreen1","darkred","gray" ))
	
	
	finalScanInFile=file.path(pathToImgFolder,"FinalScan.ini")
	
	blockPositions = parseFinalScan(finalScanInFile)
	widthO=as.numeric(as.character(blockPositions[1,2]))
	heightO=as.numeric(as.character(blockPositions[1,3]))
	centerOX=widthO/2
	centerOY=heightO/2
	widthSmallI=dim(smallImage)[1]
	heightSmallI=dim(smallImage)[2]
#draw the centers of the different clusters in the thumbnail
	for(i in 1:as.numeric(numSlides)){
		centerXSmall=centers[i,1]
		centerYSmall=round(heightSmallI/2)
		if(centerXSmall-10 >0 & centerXSmall+10<dim(smallImage)[1] & centerYSmall-10 >0 & centerYSmall+10<dim(smallImage)[2]){
			smallImage[(centerXSmall-10):(centerXSmall+10),(centerYSmall-10):(centerYSmall+10),1]=sliceColors[1,i]/255
			smallImage[(centerXSmall-10):(centerXSmall+10),(centerYSmall-10):(centerYSmall+10),2]=sliceColors[2,i]/255
			smallImage[(centerXSmall-10):(centerXSmall+10),(centerYSmall-10):(centerYSmall+10),3]=sliceColors[3,i]/255
		}
	}
	
	centers=centers[sort.list(centers),]
	
	
	
	
#assign the subimages to the sections
	blockSlice=data.frame( 2:dim(blockPositions)[1], 2:dim(blockPositions)[1], 2:dim(blockPositions)[1], 2:dim(blockPositions)[1],2:dim(blockPositions)[1], 2:dim(blockPositions)[1],stringsAsFactors=FALSE)
	
	
	for(i in 2:dim(blockPositions)[1]){
		name=as.character(blockPositions[i,1])
		x=as.numeric(as.character(blockPositions[i,2]))
		y=as.numeric(as.character(blockPositions[i,3]))
		
		cBlockXP=x/4
		cBlockYP=y/4
		posBlockX=abs(centerOX-cBlockXP)
		posBlockY=abs(centerOY-cBlockYP)
		
		posBlockXSmall=posBlockX/(widthO/widthSmallI)
		posBlockYSmall=posBlockY/(heightO/heightSmallI)
		dMin=sqrt(widthSmallI^2+heightSmallI^2)
		slice=1
		for (k in 1:as.numeric(numSlides)){
			m=rbind(c(posBlockXSmall),centers[k])
			d=dist(m)
			if(d<dMin){
				slice=k
				dMin=d
			}
		}
		
		blockSlice[(i-1),]=c(name,slice,posBlockXSmall,posBlockYSmall,x,y)
	}
#label the positions of the subimages ## not working for mac, so check for system first!!
	if(length(grep("mac",.Platform$pkgType))==0){
	    font=drawfont(weight=600,size=fontSize)
		smallImage=drawtext(smallImage,cbind(as.numeric(blockSlice[,3]),as.numeric(blockSlice[,4])),blockSlice[,1],font=font,col="black")
	}
	writeImage(smallImage,file.path(pathToOutputFolder,"smallImage.jpg"))
	colnames(blockSlice)=c("block","slice")
	sizeO=c(widthO,heightO)
	l=list(blockSlice,sizeO,smallImage)
}
