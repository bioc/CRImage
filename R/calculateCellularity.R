calculateCellularity <-
function(filename="",image=NA,classifier=NULL,cancerIdentifier=NA,KS=FALSE,maxShape=NA,minShape=NA,failureRegion=NA,colors=c(),threshold="otsu",classesToExclude=c(),numWindows=2,classifyStructures=FALSE,pixelClassifier=NA,ksToExclude=c(),densityToExclude=c(),numDensityWindows=4){
    
	if(is.null(classifier)){
		stop("Specify a classifier. No cellularity calculation applied")
	}
	if(filename !=""){
		imageData=try(segmentImage(filename=filename,image=NA,maxShape=maxShape,minShape=minShape,failureRegion=failureRegion,threshold=threshold,window=numWindows,classifyStructures=classifyStructures,pixelClassifier=pixelClassifier))
	}else{
		imageData=try(segmentImage(filename="",image=image,maxShape=maxShape,minShape=minShape,failureRegion=failureRegion,threshold=threshold,window=numWindows,classifyStructures=classifyStructures,pixelClassifier=pixelClassifier))
	}
	img=imageData[[1]]
	imgW=imageData[[2]]
	cellData=imageData[[3]]
	indexWhitePixel=imageData[[4]]
	if(classifyStructures==TRUE){
		structures=imageData[[7]]
	}
	classes=try(classifyCells(classifier=classifier,filename="",image=img,segmentedImage=imgW,featuresObjects=cellData,paint=TRUE,KS=KS,structures=structures,cancerIdentifier=cancerIdentifier,colors=colors,classifyStructures=classifyStructures))
#is this true? must have a look
	classValues=classifier$levels
	classLabels=classes[[1]]
	classLabels=classLabels[,"classCell"]
	cellularity=try(determineCellularity(classes=classLabels,classifiedCells=cellData,dimImg=dim(imgW),img=img,imgW=imgW,indexWhitePixel=indexWhitePixel,cancerIdentifier=cancerIdentifier,classValues=classValues,densityToExclude=densityToExclude,numDensityWindows=numDensityWindows))
	list(cellularity[[1]],cellularity[[2]],classes[[2]])
}
