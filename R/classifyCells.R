classifyCells <-
function(classifier,filename="",image=NA,segmentedImage=NA,featuresObjects=NA,paint=TRUE,KS=FALSE,cancerIdentifier=NA,maxShape=NA,minShape=NA,failureRegion=NA,colors=c(),classesToExclude=c(),threshold="otsu",numWindows=2,structures=NA,classifyStructures=FALSE,pixelClassifier=NA,ksToExclude=c()){
	
	message("classify cells")
	if(filename != ""){
		imageData=try(segmentImage(filename=filename,maxShape=maxShape,minShape=minShape,failureRegion=NA,threshold=threshold,numWindows=numWindows,classifyStructures=classifyStructures,pixelClassifier=pixelClassifier))
		featuresObjects=imageData[[3]]
		image=imageData[[1]]
		segmentedImage=imageData[[2]]
		if(classifyStructures==TRUE){
			structures=imageData[[7]]
		}
	}
	indexCells=featuresObjects[,c("index")]
	
	cellCoordinates=featuresObjects[,c("m.cx","m.cy")]
#if(length(classifier$x.scale[[1]])==94){
#			message("classify without topological features")
#			index=NULL
#			densityValues=NULL
#			sizeCytoplasma=NULL
#			classCell=NULL
#			g.x=NULL
#			g.y=NULL
#			g.edge=NULL
#			featuresObjectsS=subset(featuresObjects, select = -c(index,densityValues,sizeCytoplasma,classCell,m.cx,m.cy,g.edge))
#			predictedClasses=try(predict(classifier,featuresObjectsS, probability = TRUE))
#	}else{
	message("classify with topological features")
	index=NULL
	m.cx=NULL
	m.cy=NULL
	g.edge=NULL
	densityValues=NULL
	m.x=NULL
	m.y=NULL
	sizeCytoplasma=NULL
	featuresToChoose=is.element(colnames(featuresObjects),names(classifier$x.scale[[1]]))
	featuresObjects=featuresObjects[,featuresToChoose]
	
	predict(classifier,featuresObjects, probability = TRUE)
	predictedClasses=try(predict(classifier,featuresObjects, probability = TRUE))
#	}
	classSVM=predictedClasses[1:dim(featuresObjects)[1]]
	classProbs=attr(predictedClasses, "prob")
	classSVM=data.frame(indexCells,classSVM,stringsAsFactors=FALSE)
	colnames(classSVM)=c("index","classCell")
	if(KS==TRUE){
		classLabels=unique(as.character(classSVM))
		if(length(classLabels) >2 && is.na(cancerIdentifier)){
			classValues=classifier$levels
			message("More than two classes, Kernel Smoother is not applied.")
		}else if(length(classLabels) >2){
			message("Apply Kernel Smoother")
			if(is.na(cancerIdentifier)){
				message("Please specify the class for which the smoothing should be applied.")
				message("Kernel Smoother is not applied.")
			}
			classSVM=kernelSmoother(predictedClasses,cellCoordinates,segmentedImage,cancerIdentifier,indexCells,classesToExclude=classesToExclude,classValues=classifier$levels,ksToExclude=ksToExclude)
		}else{
			message("Apply Kernel Smoother with two classes")
			
			if(is.na(cancerIdentifier)){
#fix cancerIdentifier to the first class value
				cancerIdentifier=classLabels[1]
			}
#fix class to label other cells
			classSVM=kernelSmoother(predictedClasses,cellCoordinates,segmentedImage,cancerIdentifier,indexCells,classesToExclude=classesToExclude,classValues=classifier$levels,ksToExclude=ksToExclude)
			colnames(classSVM)=c("index","classCell")
		}
	}else{
		classValues=classifier$levels
	}
	
	classValues=classifier$levels
	if(paint==TRUE){
		
		paintedCells=paintCells(segmentedImage,image,as.character(classSVM[,2]),classSVM[,1],classValues,colors=colors)
		
#label all cells in structures with 1 else 0
		if(classifyStructures==TRUE){
			if(is.na(structures[1])==FALSE){
				message("Use hierarchical classification.")
				classSVMStructures=classifyStructures(structures,as.character(classSVM[,2]),classValues,classesToExclude[1],cancerIdentifier)
				classSVM=data.frame(classSVM[,1],classSVMStructures,stringsAsFactors=FALSE)
				colnames(classSVM)=c("index","classCell")
				paintedClassifyStructures=paintCells(segmentedImage,image,as.character(classSVM[,2]),classSVM[,1],classifier$levels,colors=colors)
				structures[structures != 0]=1
				paintedStructureCells=paintCells(segmentedImage,image,as.character(structures),classSVM[,1],unique(as.character(structures)),colors=colors)
				list(classSVM,paintedCells,paintedClassifyStructures,paintedStructureCells,classProbs)
			}else{
				message("No structures defined. No hierarchical classification applied.")
			}
		}else{
			list(classSVM,paintedCells,classProbs)
		}
	}else{
		list(classSVM,classProbs)
	}
	
}
