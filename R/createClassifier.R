createClassifier <-
function(trainingData,cross=FALSE,topo=TRUE){
	trainingData = trainingData[is.na(trainingData$classCell)==FALSE,]
	classes=trainingData$classCell
	if(topo==FALSE){
		index=NULL
		densityValues=NULL
		sizeCytoplasma=NULL
		classCell=NULL
		g.x=NULL
		g.y=NULL
		g.edge=NULL
		trainingData=subset(trainingData, select = -c(index,densityValues,sizeCytoplasma,classCell,g.x,g.y,g.edge))
	}else{
	index=NULL
	class=NULL
	g.x=NULL
	g.y=NULL
	g.edge=NULL
	trainingData=subset(trainingData, select = -c(index,classCell,g.x,g.y,g.edge))
	} 
	model = svm(trainingData,as.character(classes),type='C',kernel='radial',probability=TRUE)
	allCrossValues=c()
	if(cross==TRUE){
		for (i in 1:10){
			crossValue = svm(trainingData,as.character(classes),type='C',kernel='radial',probability=TRUE,cross=10)
			allCrossValues=c(allCrossValues,crossValue$a)
		}
		l=list(model,mean(allCrossValues))
	}else{
		l=list(model)
	}
	
}

