kernelSmoother <-
function(predictedClasses,cellCoordinates,segmentedImage,cancerIdentifier,indexCells,classesToExclude=c(),classValues,ksToExclude=c()){
	message("Kernel Smoother")
	prob=attr(predictedClasses, "prob")
#convert all class values to character
	predictedClasses=as.character(predictedClasses)
	classesToExclude=as.character(classesToExclude)
	ksToExclude=as.character(ksToExclude)
	classValues=as.character(classValues)
	oldClasses=predictedClasses
#allIndicesCancer=list()
#	allIndicesOther=list()
#	for (classV in 1:length(classValues)){
#		if(classValues[classV] != classOther && classValues[classV] != cancerIdentifier){
#			allIndicesOther[[classV]]=which(prob[,classOther]>=prob[,classValues[classV]])
#			allIndicesCancer[[classV]]=which(prob[,cancerIdentifier]>=prob[,classValues[classV]])
#		}
	
#	}
#	intersectCancer=allIndicesCancer[[1]]
#	for (elCL in 2:length(allIndicesCancer)){
#		if(length(allIndicesCancer[[elCL]])>0){
#			intersectCancer=intersect(intersectCancer,allIndicesCancer[[elCL]])
#		}
#	}
#	intersectOther=allIndicesOther[[1]]
#	for (elCO in 2:length(allIndicesCancer)){
#		if(length(allIndicesOther[[elCO]])>0){
#			intersectOther=intersect(intersectOther,allIndicesOther[[elCO]])
#		}
#	}
#	indicesToProcess=intersect(intersectCancer,intersectOther)
#all possible indices 
#	indicesToExclude=1:length(predictedClasses)
#	indicesToExclude=setdiff(indicesToExclude,indicesToProcess)
#these cells will get their old class after smoothing

	
	
	probClasses=prob
	cellCoordinatesN=cellCoordinates
	cellCoordinatesN=cbind(indexCells,cellCoordinates)
	cellCoordinatesN=cbind(cellCoordinatesN,1:dim(cellCoordinates)[1])
	
	
	cellCoordinatesN_temp=cellCoordinatesN
	cellCoordinatesN_temp[,4]=1:dim(cellCoordinatesN_temp)[1]
	
	
	imgWidth=dim(segmentedImage)[1]
	imgHeight=dim(segmentedImage)[2]
	xs=imgWidth/2
	xPos=imgWidth/2
	ys=imgHeight/2
	yPos=imgHeight/2
	allNewClasses=data.frame(stringsAsFactors=FALSE)
	windowCounter=0
	for (iW in 1:2){
		xPos=xs
		for (jW in 1:2){
			if(length(cellCoordinatesN_temp)>0){
				
				actualCoordinates=subset(cellCoordinatesN_temp,cellCoordinatesN_temp[,2]<=xPos & cellCoordinatesN_temp[,3]<=yPos)
				indexActualCoordinates=which(cellCoordinatesN_temp[,2]<=xPos & cellCoordinatesN_temp[,3]<=yPos)
				cellCoordinatesN_temp=subset(cellCoordinatesN_temp, !(cellCoordinatesN_temp[,4] %in% indexActualCoordinates))
				cellCoordinatesN_temp[,4]=1:dim(cellCoordinatesN_temp)[1]
				if(dim(actualCoordinates)[1]>0){
					windowCounter=windowCounter+1
					distMatrix=as.matrix(dist(actualCoordinates[,2:3]))
					fClasses=vector("list",length(classValues))
					for (i in 1:dim(actualCoordinates)[1]){
						t=distMatrix[i,]/50
						indexZero=which(t>1)
						k=(1-t^3)^3
						k[indexZero]=0
						for(classV in 1:length(classValues)){
							probClass=probClasses[,classValues[classV]]
							probClass=probClass[actualCoordinates[,"indexCells"]]
							fClasses[[classV]]=c(fClasses[[classV]],sum(k*probClass)/sum(k))
						}
					}
					index=actualCoordinates[,"indexCells"]
					smoothedClasses=fClasses
					newClassesKS=data.frame(index,stringsAsFactors=FALSE)
					for( classV in 1:length(classValues)){
						newClassesKS=cbind(newClassesKS,smoothedClasses[[classV]])
					}
					if(windowCounter==1){
						allNewClasses=newClassesKS
					}else{
						allNewClasses=rbind(allNewClasses,newClassesKS)
					}
				}
			}
			xPos=xPos+xs
		}
		yPos=yPos+ys
	}
	allNewClasses=as.data.frame(allNewClasses,stringsAsFactors=FALSE)
	
	allNewClassesN=allNewClasses
	allNewClassesN$max=apply(allNewClassesN[,2:dim(allNewClasses)[2]], 1, max)
	allNewClassesN$maxInd=apply(allNewClassesN[,2:dim(allNewClasses)[2]], 1, which.max)
	
	allNewClassesN$smoothedClass=classValues[allNewClassesN$maxInd]
	allNewClassesNKS=data.frame(allNewClassesN$index,allNewClassesN$smoothedClass,stringsAsFactors=FALSE)
	
	realIndex=data.frame(indexCells,oldClasses,stringsAsFactors=FALSE)
	colnames(realIndex)=c("index","oldClass")
	colnames(allNewClassesNKS)=c("index","classValue")
	
	
	
	classesRealIndices=merge(realIndex,allNewClassesNKS,all.x=TRUE,all.y=TRUE,stringsAsFactors=FALSE)
	
#classesRealIndices=cbind(classesRealIndices,classesRealIndices$classValues,stringsAsFactors=FALSE)
#classes which are smooted and whose value is <0.5, they keep their old label
	
	nonProbKS=which(allNewClassesN$max<0.5)
	for(classV in  classesToExclude){
		classesToReplace=which(classesRealIndices$classValue == classV)
		classesRealIndices$classValue[classesToReplace]=classesRealIndices$oldClass[classesToReplace]
#print(classesRealIndices$oldClass[classesToReplace])
	}
	indicesToReplace=c()
	if(length(ksToExclude)>0){
		#indicesToExclude=c()
		indicesToKeep=c()
		for (classE in ksToExclude){
			indicesToReplace=c(indicesToReplace,which(classesRealIndices$oldClass==classE))
		}
	}
	classesRealIndices$classValue[indicesToReplace]=classesRealIndices$oldClass[indicesToReplace]
	#print(classesRealIndices$oldClass[nonProbKS])
	classesRealIndices$classValue[nonProbKS]=classesRealIndices$oldClass[nonProbKS]
	nanIndices=which(is.na(classesRealIndices$classValue))
	#if a class was not in the window of the Kernel Smoother, replace its label by the old class
	classesRealIndices$classValue[nanIndices]=classesRealIndices$oldClass[nanIndices]
#indName=which(names(classesRealIndices)=="classesRealIndices$classValues")
#	names(classesRealIndices)[indName]="mergedClass"
#	indicesCancer=which(classesRealIndices$oldClasses == cancerIdentifier)
#	indNotCancer=classesRealIndices$mergedClass != toString(cancerIdentifier)
#	indNotCancer[is.na(indNotCancer)]=TRUE
#	classesRealIndices$mergedClass[indNotCancer]=as.character(classesRealIndices$oldClasses[indNotCancer])
#	classesRealIndices$mergedClass[indicesCancer]=as.character(classesRealIndices$classValues[indicesCancer])
#if(length(classesToExclude)>0){
#classesRealIndices[indicesToExclude,]=predictedClasses[indicesToExclude]
#}
#	allClasses=classesRealIndices[,"mergedClass"]
#	classSVM=allClasses
	classSVM=classesRealIndices[,c(1,3)]
	classSVM
}

