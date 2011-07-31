classificationAperio <-
function(fileLocation,filename,pathToOutputFolderImgDir,classifier,pathToOutputFolderImgDirFiles,pathToOutputFolderImgDirImages,blockSlice,sliceColors,sizeO,index,cancerIdentifier,maxShape=800,minShape=40,failureRegion=2000,KS=TRUE,colors=c(),classesToExclude=c(),threshold="otsu",numWindows=2,classifyStructures=FALSE,pathToOutputFolderImgDirStructures=NA,pathToOutputFolderImgDirCells,ksToExclude=c(),pixelClassifier=NA,densityToExclude=c(),numDensityWindows=32,plotCellTypeDensity=TRUE){
				widthO=as.numeric(as.character(sizeO[1]))
				heightO=as.numeric(as.character(sizeO[2]))
				pathToFile=file.path(fileLocation,filename)
				# segmentation of the image
				imageData=try(segmentImage(filename=pathToFile,maxShape=maxShape,minShape=minShape,failureRegion=NA,threshold=threshold,numWindows=numWindows,classifyStructures=classifyStructures,pixelClassifier=pixelClassifier))
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
				#name of the image
				nameFile=strsplit(filename,"\\.")[[1]][1]
				
				#is the image a failure image(no classification)
				message(nameFile)
				#section of the image
				slice=blockSlice[as.character(blockSlice$block)==nameFile,2]
				if(classify==TRUE){
					if(classifyStructures==TRUE){
						#Structure Image
						structureImages=imageData[[6]]
						structures=imageData[[7]]
						writeImage(structureImages[[1]],file.path(pathToOutputFolderImgDirStructures,paste(nameFile,"_img1.jpg",sep="")))
						#writeImage(structureImages[[2]],file.path(pathToOutputFolderImgDirStructures,paste(nameFile,"_img2.jpg",sep="")))
						#writeImage(structureImages[[3]],file.path(pathToOutputFolderImgDirStructures,paste(nameFile,"_img3.jpg",sep="")))
						#writeImage(structureImages[[4]],file.path(pathToOutputFolderImgDirStructures,paste(nameFile,"_img4.jpg",sep="")))
						#writeImage(structureImages[[5]],file.path(pathToOutputFolderImgDirStructures,paste(nameFile,"_img5_n.jpg",sep="")))
						#writeImage(structureImages[[6]],file.path(pathToOutputFolderImgDirStructures,paste(nameFile,"_img6_n.jpg",sep="")))
#print(structureImages[[5]])	
#						plot(structureImages[[5]])
					}else{
						structures=NA
					}
					#classify the cells
					classValues=classifyCells(classifier,image=img,segmentedImage=imgW,featuresObjects=cellData,paint=TRUE,KS=KS,cancerIdentifier=cancerIdentifier,colors=colors,classesToExclude=classesToExclude,structures=structures,classifyStructures=classifyStructures,ksToExclude=ksToExclude)
					classes=classValues[[1]]

					writeImage(classValues[[2]],file.path(pathToOutputFolderImgDir,file.path(paste("section",slice,sep="_"),filename)))
					if(classifyStructures==TRUE){
						writeImage(classValues[[3]],file.path(pathToOutputFolderImgDir,file.path(paste("section",slice,sep="_"),paste(nameFile,"_structure.jpg",sep=""))))
						classProbs=classValues[[4]]
					}else{
						classProbs=classValues[[3]]
					}
					cellData=merge(classes,cellData)
					paintedCellsFinal=paintCells(imgW,img,as.character(cellData[,"classCell"]),cellData[,"index"],classifier$levels,colors=colors)
					#writeImage(paintedCellsFinal,file.path(pathToOutputFolderImgDir,file.path(paste("section",slice,sep="_"),paste(nameFile,"_Final.jpg"))))

					#calculation of cellularity, drawing of the heatmap
					cellValues=determineCellularity(cellData[,"classCell"],cellData,dim(imgW),img,imgW,indexWhitePixel,cancerIdentifier,classifier$levels,densityToExclude=densityToExclude,numDensityWindows=numDensityWindows,plotCellTypeDensity=plotCellTypeDensity)
					cellsSubImages=cellValues[[1]]
					cellsDensityImage=cellValues[[2]]
					
					cellValueTable=data.frame(t(cellsSubImages),stringsAsFactors=FALSE)
					colnames(cellValueTable)=names(cellsSubImages)
					write.table(cellValueTable,file.path(pathToOutputFolderImgDirFiles,file.path(paste("section",slice,sep="_"),nameFile)),sep="\t",row.names = FALSE)
					cellDataProbs=data.frame(cellData,classProbs,stringsAsFactors=FALSE)	
					write.table(cellDataProbs,file.path(pathToOutputFolderImgDirCells,paste(nameFile,".txt",sep="")),sep="\t",row.names=FALSE,quote = FALSE)
					message("Cell values written")	
					message("Write tumour density heatmap")
					writeDensityImage(pathToOutputFolderImgDir,cellsDensityImage,sizeO,blockSlice,"smallDensityImage.jpg",nameFile)
					
					if(plotCellTypeDensity==TRUE){
						cellsTypeImage=cellValues[[3]]
						writeDensityImage(pathToOutputFolderImgDir,cellsTypeImage,sizeO,blockSlice,"cellTypeImage.jpg",nameFile)
					}
				###########################################################################################################################			

				
					cellData
				}else{
					
###########################################################################################################################			
					
				classifiedCells=data.frame(stringsAsFactors=FALSE)
				}
}

