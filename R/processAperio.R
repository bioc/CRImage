processAperio <-
function(classifier=classifier,inputFolder=inputFolder,outputFolder=outputFolder,identifier=identifier,numSlides=numSlides,cancerIdentifier=cancerIdentifier,maxShape=NA,minShape=NA,failureRegion=NA,slideToProcess=NA,KS=TRUE){
	options(stringsAsFactors = FALSE)
	pathToFolder=inputFolder
	pathToOutputFolder=outputFolder
	pathToOutputFolderTempFiles=file.path(pathToOutputFolder,"tempFiles")
	dir.create(pathToOutputFolderTempFiles)
	
	dir.create(pathToOutputFolder)
	pathToOutputFolderImgDir=file.path(pathToOutputFolder,"classifiedImage")
	pathToOutputFolderImgDirFiles=file.path(pathToOutputFolder,"Files")
	pathToOutputFolderImgDirCellDensity=file.path(pathToOutputFolder,"tumourDensity")
	dir.create(pathToOutputFolderImgDir)
	dir.create(pathToOutputFolderImgDirFiles)
	dir.create(pathToOutputFolderImgDirCellDensity)
	#finds the right section for every subimage
	sliceSizeList=try(findSlices(inputFolder,pathToOutputFolder,numSlides))
	blockSlice=sliceSizeList[[1]]
	sizeO=sliceSizeList[[2]]
	numberSlices=length(unique(blockSlice[,2]))
	sliceFolder=c()
	for (i in 1:numberSlices){
		dir.create(file.path(pathToOutputFolderImgDirFiles,paste("section",i,sep="_")))
	}
	sliceColors=col2rgb(c("red","blue","green","yellow","orange"))
	dir.create(pathToOutputFolderImgDirCellDensity)
	filenames=list.files(path =pathToFolder ,pattern=identifier)
	allCells=foreach( i = 1:length(filenames), .combine=rbind) %do% {
	nameFile=strsplit(filenames[i],"\\.")[[1]][1]
		if(is.na(slideToProcess)){
			classificationError=try(classificationAperio(pathToFolder,filenames[i],pathToOutputFolderImgDir,classifier,pathToOutputFolderImgDirFiles,pathToOutputFolderImgDir,pathToOutputFolderImgDirCellDensity,blockSlice,sliceColors,sizeO,i,cancerIdentifier,KS=KS))
		}else{
			slice=blockSlice[as.character(blockSlice$block)==nameFile,2]
			if(slice == slideToProcess){
				classificationError=try(classificationAperio(pathToFolder,filenames[i],pathToOutputFolderImgDir,classifier,pathToOutputFolderImgDirFiles,pathToOutputFolderImgDir,pathToOutputFolderImgDirCellDensity,blockSlice,sliceColors,sizeO,i,cancerIdentifier,KS=KS))
			}else{
				message(paste("Subimage",paste(nameFile, "is not processed.")))
			}
				
		}
						}
		write.table(allCells,file.path(pathToOutputFolderImgDir,paste("result",".txt",sep="")),col.names=TRUE,sep="\t",row.names=FALSE)
}

