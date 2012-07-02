determineCellularity <-
function(classes,classifiedCells,dimImg,img,imgW,indexWhitePixel,cancerIdentifier,classValues,densityToExclude=c(),numDensityWindows=32,plotCellTypeDensity=TRUE){
	message("calculate Cellularity.")
	wholeCellDensityImage=img
	if(plotCellTypeDensity==TRUE){
		wholeCellTypeImage=img
	}
	imgTC=imgW
	imgW[img[,,1]==2]=-1
	rm(img)
	imgW[indexWhitePixel]=-1
	predictedClassesN=c(as.character(classes),cancerIdentifier)
	imgTC[imgTC==0]=length(predictedClassesN)
	a=array(predictedClassesN[imgTC] != cancerIdentifier,dim(imgW))
	imgTC[a]=0
	imgTC[imgTC==length(predictedClassesN)]=0
	
	classifiedCells[,"g.x"]=as.numeric(as.character(classifiedCells[,"g.x"]))
	classifiedCells[,"g.y"]=as.numeric(as.character(classifiedCells[,"g.y"]))
	numWindows=numDensityWindows
	xStepSize=dimImg[1]/numWindows
	yStepSize=dimImg[2]/numWindows
	xl=1
	xr=xStepSize
	yo=1
	yu=yStepSize
#hColors=col2rgb(heat.colors(10))
#hColors=c(1,25,50,75,100,125,150,175,200,255)
#hColors=hColors[,dim(hColors)[2]:1]
#	rColors=hColors[1,]/256
#	gColors=hColors[2,]/256
#	bColors=hColors[3,]/256
	
	
	s=seq(1,255,20)
	cellTypeColors=s/255
	
	rgb.palette = colorRampPalette(c("white", "lightgreen", "green","darkgreen"),space = "rgb")
	
	HColors=col2rgb(rgb.palette(20))
	rColors=HColors[1,]/255
	gColors=HColors[2,]/255
	bColors=HColors[3,]/255
	
	
	cellularity=length(classes[classes==cancerIdentifier])/length(imgW[imgW != -1])
	numTumorPixel=length(imgW[imgW != -1])
	cellularityValues=c()
	numClassCells=c()
	cancerCells=0
	numRealCells=c()
	for (classValue in classValues){
		excludeThisClass=FALSE
		for(classE in densityToExclude){
			if (classE == classValue){
				excludeThisClass=TRUE
			}
		}
		if (excludeThisClass==FALSE){
			numRealCells=c(numRealCells,length(classes[classes==classValue]))
		}
		numClassCells=c(numClassCells,length(classes[classes==classValue]))
		if(classValue == cancerIdentifier){
			cancerCells=length(classes[classes==classValue])
		}
	}
	
	indexValues=c(classValues,"cellularity","ratioTumourCells","numTumourPixel")
	
	
	numAllCells=sum(numClassCells)
	ratioCancerCells=cancerCells/sum(numRealCells)
	cellularityValues=c(cellularityValues,numClassCells,cellularity,ratioCancerCells,numTumorPixel)
	names(cellularityValues)=indexValues
	cancerCells=c()
	for (i in 1:numWindows){
		xl=1
		xr=xStepSize
		for (j in 1:numWindows){
			
			imgSub1=imgW[xl:xr,yo:yu]
			numPositivePixels=length(imgSub1[imgSub1 != -1])						
			cellsWindow=c()
			cancerCells=c()
			otherCellsWindow=c()
			for (classValue in classValues){
				cells=subset(classifiedCells,classifiedCells[,"g.x"]<=xr & classifiedCells[,"g.x"]>=xl & classifiedCells[,"g.y"]<=yu & classifiedCells[,"g.y"]>=yo & classes==classValue)
				excludeThisClass=FALSE
				for(classE in densityToExclude){
					if (classE == classValue){
						excludeThisClass=TRUE
					}
				}
				if (excludeThisClass==FALSE){
					otherCellsWindow=c(otherCellsWindow,dim(cells)[1])
				}
				
				if(classValue == cancerIdentifier){
					cancerCells=cells
				}
				
			}
			if(length(cancerCells)==0){
				ratioCancerCellPixel=0
			}else{
#ratioCancerCellPixel=dim(cancerCells)[1]/numPositivePixels
				ratioCancerCellPixel=dim(cancerCells)[1]/sum(otherCellsWindow)
				
			}
			cellRatios=rep(0,length(classValues))
			if(length(otherCellsWindow)<=3){
				for (cellV in 1:length(otherCellsWindow)){
					if(!is.na(sum(otherCellsWindow))){
						if(sum(otherCellsWindow)!=0){
							cellRatios[cellV]=otherCellsWindow[cellV]/sum(otherCellsWindow)
						}
					}
				}
			}
			if(is.na(ratioCancerCellPixel)){
				ratioCancerCellPixel=0
			}
#colorRatioCancerCellPixel=(ratioCancerCellPixel*500)*50+1
			colorRatioCancerCellPixel=ceiling(ratioCancerCellPixel*length(rColors))
			if(colorRatioCancerCellPixel==0){
				colorRatioCancerCellPixel=1
			}
			if(colorRatioCancerCellPixel>length(rColors)){
				colorRatioCancerCellPixel=length(rColors)
			}
			if(is.na(ratioCancerCellPixel)){
				wholeCellDensityImage[xl:xr,yo:yu,1]=1
				wholeCellDensityImage[xl:xr,yo:yu,2]=1
				wholeCellDensityImage[xl:xr,yo:yu,3]=1
			}else{
				
#print(cellRatios)
#				print(ceiling(cellRatios[3]*length(cellTypeColors)))
#				print(ceiling(cellRatios[1]*length(cellTypeColors)))
#				print(ceiling(cellRatios[2]*length(cellTypeColors)))
				
				
				wholeCellDensityImage[xl:xr,yo:yu,1]=rColors[colorRatioCancerCellPixel]
				wholeCellDensityImage[xl:xr,yo:yu,2]=gColors[colorRatioCancerCellPixel]
				wholeCellDensityImage[xl:xr,yo:yu,3]=bColors[colorRatioCancerCellPixel]
				
				if(plotCellTypeDensity==TRUE && length(cellRatios)>=3){
					if(cellRatios[1]==0 && cellRatios[2]==0 && cellRatios[3]==0){
						cellRatios[1]=1
						cellRatios[2]=1
						cellRatios[3]=1
					}
					if(cellRatios[1]==0){
						cellRatios[1]=0.001
					}
					if(cellRatios[2]==0){
						cellRatios[2]=0.001
					}
					if(cellRatios[3]==0){
						cellRatios[3]=0.001
					}
					wholeCellTypeImage[xl:xr,yo:yu,1]=cellTypeColors[ceiling(cellRatios[3]*length(cellTypeColors))]
					wholeCellTypeImage[xl:xr,yo:yu,2]=cellTypeColors[ceiling(cellRatios[1]*length(cellTypeColors))]
					wholeCellTypeImage[xl:xr,yo:yu,3]=cellTypeColors[ceiling(cellRatios[2]*length(cellTypeColors))]
				}
			}
			xl=xl+xStepSize
			xr=xr+xStepSize
			
		}
		yo=yo+yStepSize
		yu=yu+yStepSize
	}
	
	wholeCellDensityImage=Image(wholeCellDensityImage)
	colorMode(wholeCellDensityImage)=Color
	
	message("density subimage created")
	if(plotCellTypeDensity==TRUE){
		wholeCellTypeImage=Image(wholeCellTypeImage)
		colorMode(wholeCellTypeImage)=Color
		l=list(cellularityValues,wholeCellDensityImage,wholeCellTypeImage)
	}else{
		l=list(cellularityValues,wholeCellDensityImage)
	}
	
}
