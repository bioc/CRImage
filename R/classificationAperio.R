classificationAperio <-
function(fileLocation,filename,pathToOutputFolderImgDir,classifier,pathToOutputFolderImgDirFiles,pathToOutputFolderImgDirImages,pathToOutputFolderImgDirCellDensity,blockSlice,sliceColors,sizeO,index,cancerIdentifier,maxShape=NA,minShape=NA,failureRegion=NA,KS=TRUE){
				widthO=as.numeric(as.character(sizeO[1]))
				heightO=as.numeric(as.character(sizeO[2]))
				pathToFile=file.path(fileLocation,filename)
				# segmentation of the image
				imageData=try(segmentImage(filename=pathToFile,maxShape=800,minShape=40,failureRegion=2000))
				#original image
				img=imageData[[1]]
				#segmented image
				imgW=imageData[[2]]
				#features of the nulcei
				cellData=imageData[[3]]
				#indices of all white pixel in the image
				indexWhitePixel=imageData[[4]]
				#classify the nuclei, or is the image white?
				classify=imageData[[5]]
				#is the image a failure image(no classification)
				nameFile=strsplit(filename,"\\.")[[1]][1]
				#section of the image
				slice=blockSlice[as.character(blockSlice$block)==nameFile,2]
				if(classify==TRUE){
			    #classify the cells
				classValues=classifyCells(classifier,image=img,segmentedImage=imgW,featuresObjects=cellData,paint=TRUE,KS=KS,cancerIdentifier=cancerIdentifier)
				classes=classValues[[1]]
				writeImage(classValues[[2]],file.path(pathToOutputFolderImgDir,filename))
				classes=as.character(classes)
				classifiedCells=cellData
				#calculation of cellularity, drawing of the heatmap
				cellValues=determineCellularity(classes,classifiedCells,dim(imgW),img,imgW,indexWhitePixel,cancerIdentifier,sort(unique(classes)))
				cellsSubImages=cellValues[[1]]
				cellsDensityImage=cellValues[[2]]
				cellValueTable=data.frame(t(cellsSubImages))
				colnames(cellValueTable)=names(cellsSubImages)
				write.table(cellValueTable,file.path(pathToOutputFolderImgDirFiles,file.path(paste("section",slice,sep="_"),nameFile)),sep="\t",row.names = FALSE)
				
				#draw the heatmap
				writeImage(cellsDensityImage,file.path(pathToOutputFolderImgDirCellDensity,paste(nameFile,".jpg",sep="")))
				classesN=classes
				classesN[which(classes==cancerIdentifier)]=1
				classesN[which(classes!=cancerIdentifier)]=2
				classifiedCells[,"classCell"]=as.numeric(classesN)
				classifiedCells
				}else{
					colorMode(img)="color"	
					nameFile=strsplit(filename,"\\.")[[1]][1]
					if(length(dim(img))>2){
						cellsDensityImage=img
					}else{
						cellsDensityImage=img
					}
					if(length(dim(img))>2){
						cellsDensityImage[,,1]=1
						cellsDensityImage[,,2]=1
						cellsDensityImage[,,3]=1
					}
					writeImage(cellsDensityImage,file.path(pathToOutputFolderImgDirCellDensity,paste(nameFile,".jpg",sep="")))
						
					writeImage(img,file.path(pathToOutputFolderImgDir,paste(nameFile,".jpg",sep="")))
					
					classifiedCells=data.frame()
				}
}

