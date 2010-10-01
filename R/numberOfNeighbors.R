numberOfNeighbors <-
function(img,cellCoordinates,allFeatures){
	xs=dim(img)[1]/2
	xPos=dim(img)[1]/2
	ys=dim(img)[2]/2
	yPos=dim(img)[2]/2
	cellCoordinatesN=cellCoordinates
	cellCoordinatesN=cbind(allFeatures[,"index"],cellCoordinates)
	cellCoordinatesN=cbind(cellCoordinatesN,1:dim(cellCoordinates)[1])
	#calculate the number of neighbors of every nuclei
	indexNeighbors=data.frame()
	for (i in 1:2){
		xPos=xs
		for (j in 1:2){
			if(length(cellCoordinatesN)>0){
				actualCoordinates=subset(cellCoordinatesN,cellCoordinatesN[,2]<=xPos & cellCoordinatesN[,3]<=yPos)
				indexActualCoordinates=which(cellCoordinatesN[,2]<=xPos & cellCoordinatesN[,3]<=yPos)
				cellCoordinatesN=subset(cellCoordinatesN, !(cellCoordinatesN[,4] %in% indexActualCoordinates))
				cellCoordinatesN[,4]=1:dim(cellCoordinatesN)[1]
				if(dim(actualCoordinates)[1]>0){
					distMatrix=as.matrix(dist(actualCoordinates[,2:3]))
					numberNeighbors=c()
					f=c()
					for (m in 1:dim(actualCoordinates)[1]){
						sortedCells=sort(distMatrix[m,])
						neighbors=length(sortedCells[sortedCells<50])
						numberNeighbors=c(numberNeighbors,neighbors)
					}
					numberNeighbors=cbind(actualCoordinates[,1],numberNeighbors)
					indexNeighbors=rbind(indexNeighbors,numberNeighbors)
				}
			}
			xPos=xPos+xs
		}
		yPos=yPos+ys
	}
	realIndex=data.frame(allFeatures[,"index"])
	colnames(realIndex)="index"
	colnames(indexNeighbors)=c("index","neighbors")
	neighborsRealIndices=merge(realIndex,indexNeighbors,all.x=TRUE,all.y=TRUE)
	numberNeighbors=neighborsRealIndices$neighbors
	numberNeighbors
}

