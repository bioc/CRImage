createClassifier <-
function(trainingData,cross=FALSE){
	#trainingData = trainingData[is.na(trainingData$classCell)==FALSE,]
	classes=trainingData$classCell
		index=NULL
		densityValues=NULL
		sizeCytoplasma=NULL
		classCell=NULL
		g.x=NULL
		g.y=NULL
		g.edge=NULL
		#trainingData=subset(trainingData, select = -c(index,densityValues,sizeCytoplasma,classCell,g.x,g.y,g.edge))
		indToDelete=which(is.element(colnames(trainingData),c("index","classCell","m.cx","m.cy","g.edge","sizeCytoplasma","densityValues")))
		trainingData=subset(trainingData,select=-c(indToDelete))
		
	model = svm(trainingData,as.factor(classes),type='C',kernel='radial',probability=TRUE)
	allCrossValues=c()
	if(cross==TRUE){
		for (i in 1:10){
			crossValue = svm(trainingData,as.factor(classes),type='C',kernel='radial',probability=TRUE,cross=10)
			allCrossValues=c(allCrossValues,crossValue$a)
		}
		l=list(model,mean(allCrossValues))
	}else{
		l=list(model)
	}
	
}
