classifyCells <-
function(classifier,filename="",image=NA,segmentedImage=NA,featuresObjects=NA,paint=TRUE,KS=FALSE,cancerIdentifier=NA,maxShape=NA,minShape=NA,failureRegion=NA){
	if(filename != ""){
		imageData=try(segmentImage(filename,maxShape=maxShape,minShape=minShape,failureRegion=failureRegion))
		featuresObjects=imageData[[3]]
		image=imageData[[1]]
		segmentedImage=imageData[[2]]
	}
	indexCells=featuresObjects[,c("index")]
	
	cellCoordinates=featuresObjects[,c("g.x","g.y")]
	if(length(classifier$x.scale[[1]])==94){
			message("classify without topological features")
			index=NULL
			densityValues=NULL
			sizeCytoplasma=NULL
			classCell=NULL
			g.x=NULL
			g.y=NULL
			g.edge=NULL
			featuresObjectsS=subset(featuresObjects, select = -c(index,densityValues,sizeCytoplasma,classCell,g.x,g.y,g.edge))
			predictedClasses=try(predict(classifier,featuresObjectsS, probability = TRUE))
	}else{
			message("classify with topological features")
			index=NULL
			classCell=NULL
			g.x=NULL
			g.y=NULL
			g.edge=NULL
			featuresObjectsS=subset(featuresObjects, select = -c(index,classCell,g.x,g.y,g.edge))
			predictedClasses=try(predict(classifier,featuresObjectsS, probability = TRUE))
	}
	classSVM=predictedClasses[1:dim(featuresObjects)[1]]
	if(KS==TRUE){
		classLabels=unique(as.character(classSVM))
		if(length(classLabels) >2 && is.na(cancerIdentifier)){
			classValues=classifier$levels
			message("More than two classes, Kernel Smoother is not applied.")
		}else if(length(classLabels) >2 && is.na(cancerIdentifier)==FALSE){
			message("Apply Kernel Smoother")
			classSVM=kernelSmoother(predictedClasses,cellCoordinates,segmentedImage,cancerIdentifier,indexCells)
			classValues=sort(unique(classSVM))
		}else{
			message("Apply Kernel Smoother")
			classSVM=kernelSmoother(predictedClasses,cellCoordinates,segmentedImage,1,indexCells)
			classValues=sort(unique(classSVM))
		}
	}else{
		classValues=classifier$levels
	}
		if(paint==TRUE){
			paintedCells=paintCells(segmentedImage,image,as.character(classSVM),indexCells,classValues)
			list(as.character(classSVM),paintedCells)
		}else{
			list(as.character(classSVM))
		}

}

