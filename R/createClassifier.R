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
			m.x=NULL
			m.y=NULL
#trainingData=subset(trainingData, select = -c(index,densityValues,sizeCytoplasma,classCell,g.x,g.y,g.edge))
			trainingData=subset(trainingData,select=-c(index,classCell,g.x,g.y,g.edge,m.x,m.y,sizeCytoplasma,densityValues))
		}else{
			index=NULL
			class=NULL
			g.x=NULL
			g.y=NULL
			g.edge=NULL
			m.x=NULL
			m.y=NULL
#trainingData=subset(trainingData, select = -c(index,classCell,g.x,g.y,g.edge))
			trainingData=subset(trainingData,select=-c(index,classCell,g.x,g.y,g.edge,m.x,m.y,sizeCytoplasma,densityValues))
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
