writeDensityImage <-
function(pathToOutputFolderImgDir,cellsDensityImage,sizeO,blockSlice,name,nameFile){
	message("Write density image")
###### insert densityImage in smallDensityImage
	smallDensityImage=readImage(file.path(pathToOutputFolderImgDir,name))	
	if(length(dim(smallDensityImage))==2){
		smallDensityImage=rgbImage(red=smallDensityImage,green=smallDensityImage,blue=smallDensityImage)
	}
#smallCellsDensityImageWidth=floor(dim(cellsDensityImage)[1]*(dim(smallDensityImage)[1]/sizeO[1]))
#					smallCellsDensityImageHeight=floor(dim(cellsDensityImage)[2]*(dim(smallDensityImage)[2]/sizeO[2]))
	smallCellsDensityImageWidth=ceiling(dim(cellsDensityImage)[1]*(dim(smallDensityImage)[1]/sizeO[1]))
	smallCellsDensityImageHeight=ceiling(dim(cellsDensityImage)[2]*(dim(smallDensityImage)[2]/sizeO[2]))
	if(smallCellsDensityImageWidth>0 && smallCellsDensityImageHeight>0){
		
		smallCellsDensityImage=resize(cellsDensityImage,smallCellsDensityImageWidth,smallCellsDensityImageHeight)
		xPosition=as.numeric(blockSlice[as.character(blockSlice$block)==nameFile,5])
		yPosition=as.numeric(blockSlice[as.character(blockSlice$block)==nameFile,6])
#sizeX=floor(smallCellsDensityImageWidth/2)
#						sizeY=floor(smallCellsDensityImageHeight/2)
		sizeX=floor(smallCellsDensityImageWidth/2)
		sizeY=floor(smallCellsDensityImageHeight/2)
		sizeXRight=sizeX
		sizeYDown=sizeY
		if(smallCellsDensityImageWidth%%2 == 0){
			sizeXRight=sizeXRight-1
		}
		if(smallCellsDensityImageHeight%%2 == 0){
			sizeYDown=sizeYDown-1
		}
		
################################
		
####### 2000 ersetzen #####
		xMovement=(2000/2)-dim(cellsDensityImage)[1]/2
		yMovement=(2000/2)-dim(cellsDensityImage)[2]/2
		centerOX=sizeO[1]/2
		centerOY=sizeO[2]/2
		cBlockXP=xPosition/4
		cBlockYP=yPosition/4
		if(cBlockXP<0){
			cBlockXP=cBlockXP+xMovement
		}else{
			cBlockXP=cBlockXP-xMovement
		}
		if(cBlockYP<0){
			cBlockYP=cBlockYP+yMovement
		}else{
			cBlockYP=cBlockYP-yMovement
		}
		
		posBlockX=abs(centerOX-cBlockXP)
		posBlockY=abs(centerOY-cBlockYP)
		
		xPosition=ceiling(posBlockX*(dim(smallDensityImage)[1]/sizeO[1]))
		yPosition=ceiling(posBlockY*(dim(smallDensityImage)[2]/sizeO[2]))

		if((xPosition-sizeX)<=0){
			sizeX=sizeX-1
			smallCellsDensityImage=smallCellsDensityImage[2:dim(smallCellsDensityImage)[1],,]
		}
		if((xPosition+sizeX)>dim(smallDensityImage)[2]){
			sizeXRight=sizeXRight-1
			smallCellsDensityImage=smallCellsDensityImage[1:(dim(smallCellsDensityImage)[1]-1),,]
		}
		if((yPosition-sizeY)<=0){
#			counterUp=1
#			while((yPosition-sizeY)<=0){
				sizeY=sizeY-1
#				counterUp=counterUp+1
#}
			
			smallCellsDensityImage=smallCellsDensityImage[,2:dim(smallCellsDensityImage)[2],]
		}
		if((yPosition+sizeYDown)>dim(smallDensityImage)[2]){
			sizeYDown=sizeYDown-1
			smallCellsDensityImage=smallCellsDensityImage[,1:(dim(smallCellsDensityImage)[2]-1),]
		}
		

		imgN=smallDensityImage[(xPosition-sizeX):(xPosition+sizeXRight),(yPosition-sizeY):(yPosition+sizeYDown),]
		
		smallDensityImage[(xPosition-sizeX):(xPosition+sizeXRight),(yPosition-sizeY):(yPosition+sizeYDown),]=smallCellsDensityImage
#draw the heatmap
		writeImage(smallDensityImage,file.path(pathToOutputFolderImgDir,name))
	}
}

