processAperio <-
function(classifier=classifier,inputFolder=inputFolder,outputFolder=outputFolder,identifier=identifier,numSlides=numSlides,cancerIdentifier=cancerIdentifier,classOther=NA,maxShape=800,minShape=40,failureRegion=2000,slideToProcess=NA,KS=TRUE,colors=c(),classesToExclude=c(),threshold="otsu",window=2,colorCorrection=FALSE,classifyStructures=FALSE,ksToExclude=c(),pixelClassifier=NA,densityToExclude=c(),numDensityWindows=32,resizeFactor=4,plotCellTypeDensity=TRUE,greyscaleImage=0,penClassifier=NULL,referenceHist=NULL,fontSize=10){
	options(stringsAsFactors = FALSE)
	pathToFolder=inputFolder
	pathToOutputFolder=outputFolder
	pathToOutputFolderTempFiles=file.path(pathToOutputFolder,"tempFiles")
	dir.create(pathToOutputFolderTempFiles)
	
	dir.create(pathToOutputFolder)
	pathToOutputFolderImgDir=file.path(pathToOutputFolder,"classifiedImage")
	
	pathToOutputFolderImgDirFiles=file.path(pathToOutputFolder,"Files")
	pathToOutputFolderImgDirStructures=file.path(pathToOutputFolder,"Structures")
	pathToOutputFolderImgDirCells=file.path(pathToOutputFolder,"Cells")
	dir.create(pathToOutputFolderImgDir)
	dir.create(pathToOutputFolderImgDirFiles)
	dir.create(pathToOutputFolderImgDirStructures)
	dir.create(pathToOutputFolderImgDirCells)
#finds the right section for every subimage
	sliceSizeList=try(findSlices(inputFolder,pathToOutputFolder,numSlides,fontSize=fontSize))
	blockSlice=sliceSizeList[[1]]
	sizeO=sliceSizeList[[2]]
	smallImage=sliceSizeList[[3]]
	
	densityImage=mat.or.vec(round(dim(smallImage)[1]*resizeFactor),round(dim(smallImage)[2]*resizeFactor))
	densityImage[,]=1
	densityImageRGB=rgbImage(red=densityImage,green=densityImage,blue=densityImage)
	densityImageRGB[,,1]=1
	densityImageRGB[,,2]=1
	densityImageRGB[,,3]=1
	writeImage(densityImageRGB,file.path(pathToOutputFolderImgDir,"smallDensityImage.jpg"))
	if(plotCellTypeDensity==TRUE){
		writeImage(densityImageRGB,file.path(pathToOutputFolderImgDir,"cellTypeImage.jpg"))
	}
	qualityTable=t(c("NA","NA","NA"))
	write.table(qualityTable,file.path(pathToOutputFolderImgDir,"subimageQuality.txt"),sep="\t",row.names = FALSE)
	numberSlices=length(unique(blockSlice[,2]))
	sliceFolder=c()
	for (i in 1:numberSlices){
		dir.create(file.path(pathToOutputFolderImgDirFiles,paste("section",i,sep="_")))
		dir.create(file.path(pathToOutputFolderImgDir,paste("section",i,sep="_")))
	}
	sliceColors=col2rgb(c("red","blue","green","yellow","orange"))
	filenames=list.files(path =pathToFolder ,pattern=identifier)
	allCells=foreach( i = 1:length(filenames), .combine=rbind) %do% {
		nameFile=strsplit(filenames[i],"\\.")[[1]][1]
		if(is.na(slideToProcess)){
			classificationError=try(classificationAperio(pathToFolder,filenames[i],pathToOutputFolderImgDir,classifier,pathToOutputFolderImgDirFiles,pathToOutputFolderImgDir,blockSlice,sliceColors,sizeO,i,cancerIdentifier,maxShape=maxShape,minShape=minShape,failureRegion=failureRegion,KS=KS,colors=colors,classesToExclude=classesToExclude,threshold=threshold,window=window,classOther=classOther,pathToOutputFolderImgDirStructures,colorCorrection=colorCorrection,classifyStructures=classifyStructures,pathToOutputFolderImgDirCells,ksToExclude=ksToExclude,pixelClassifier=pixelClassifier,densityToExclude=densityToExclude,numDensityWindows=numDensityWindows,plotCellTypeDensity=plotCellTypeDensity,greyscaleImage=greyscaleImage,penClassifier=penClassifier,referenceHist=referenceHist))
		}else{
			slice=blockSlice[as.character(blockSlice$block)==nameFile,2]
			if(slice == slideToProcess){
				classificationError=try(classificationAperio(pathToFolder,filenames[i],pathToOutputFolderImgDir,classifier,pathToOutputFolderImgDirFiles,pathToOutputFolderImgDir,blockSlice,sliceColors,sizeO,i,cancerIdentifier,maxShape=maxShape,minShape=minShape,failureRegion=failureRegion,KS=KS,colors=colors,classesToExclude=classesToExclude,threshold=threshold,window=window,classOther=classOther,pathToOutputFolderImgDirStructures,colorCorrection=colorCorrection,classifyStructures=classifyStructures,pathToOutputFolderImgDirCells,ksToExclude=ksToExclude,pixelClassifier=pixelClassifier,densityToExclude=densityToExclude,numDensityWindows=numDensityWindows,plotCellTypeDensity=plotCellTypeDensity,greyscaleImage=greyscaleImage,penClassifier=penClassifier,referenceHist=referenceHist))
			}else{
				message(paste("Subimage",paste(nameFile, "is not processed.")))
			}
			
		}
	}
	write.table(allCells,file.path(pathToOutputFolderImgDir,paste("result",".txt",sep="")),col.names=TRUE,sep="\t",row.names=FALSE)
}
