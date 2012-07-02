classifyCells <-
function(classifier,filename="",image=NA,segmentedImage=NA,featuresObjects=NA,paint=TRUE,KS=FALSE,cancerIdentifier=NA,classOther=NA,maxShape=NA,minShape=NA,failureRegion=NA,colors=c(),classesToExclude=c(),threshold="otsu",window=2,structures=NA,classifyStructures=FALSE,pixelClassifier=NA,ksToExclude=c()){
	message("classify cells")
	if(filename != ""){
		imageData=try(segmentImage(filename,maxShape=maxShape,minShape=minShape,failureRegion=failureRegion,threshold=threshold,window=window,colorCorrection=FALSE))
		featuresObjects=imageData[[3]]
		image=imageData[[1]]
		segmentedImage=imageData[[2]]
	}
	indexCells=featuresObjects[,c("index")]
	
	cellCoordinates=featuresObjects[,c("g.x","g.y")]
#if(length(classifier$x.scale[[1]])==94){
#			message("classify without topological features")
#			index=NULL
#			densityValues=NULL
#			sizeCytoplasma=NULL
#			classCell=NULL
#			g.x=NULL
#			g.y=NULL
#			g.edge=NULL
#			featuresObjectsS=subset(featuresObjects, select = -c(index,densityValues,sizeCytoplasma,classCell,g.x,g.y,g.edge,m.x,m.y))
#			predictedClasses=try(predict(classifier,featuresObjectsS, probability = TRUE))
#	}else{
	message("classify with topological features")
	index=NULL
	g.x=NULL
	g.y=NULL
	g.edge=NULL
	#featuresObjectsS=subset(featuresObjects, select = -c(densityValues,index,g.x,g.y,g.edge,m.x,m.y,sizeCytoplasma))
	featuresObjectsS=subset(featuresObjects, select = c(names(classifier$x.scale[[1]])))
	predictedClasses=try(predict(classifier,featuresObjectsS, probability = TRUE))
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
			else if(is.na(classOther)){
				message("Please specify the class with which changed classes should be labelled.")
				message("Kernel Smoother is not applied.")
			}else{
				classSVM=kernelSmoother(predictedClasses,cellCoordinates,segmentedImage,cancerIdentifier,indexCells,classesToExclude=classesToExclude,classValues=classifier$levels,classOther=classOther,ksToExclude=ksToExclude)
			}
		}else{
			message("Apply Kernel Smoother with two classes")
			
			if(is.na(cancerIdentifier)){
#fix cancerIdentifier to the first class value
				cancerIdentifier=classLabels[1]
			}
#fix class to label other cells
			if(is.na(classOther)){
				classOther=classLabels[2]
			}
			classSVM=kernelSmoother(predictedClasses,cellCoordinates,segmentedImage,cancerIdentifier,indexCells,classesToExclude=classesToExclude,classValues=classifier$levels,classOther=classOther,ksToExclude=ksToExclude)
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
			classSVMStructures=classifyStructures(structures,as.character(classSVM[,2]),classValues,classesToExclude[1],cancerIdentifier,classOther)
			classSVM=data.frame(classSVM[,1],classSVMStructures,stringsAsFactors=FALSE)
			colnames(classSVM)=c("index","classCell")
			paintedClassifyStructures=paintCells(segmentedImage,image,as.character(classSVM[,2]),classSVM[,1],classifier$levels,colors=colors)
			structures[structures != 0]=1
			paintedStructureCells=paintCells(segmentedImage,image,as.character(structures),classSVM[,1],unique(as.character(structures)),colors=colors)
			list(classSVM,paintedCells,paintedStructureCells,paintedClassifyStructures,classProbs)
		}else{
			list(classSVM,paintedCells,classProbs)
		}
	}else{
		list(classSVM,classProbs)
	}
	
}
