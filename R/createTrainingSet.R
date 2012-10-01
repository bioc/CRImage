createTrainingSet <-
function(filename="",image=NA,maxShape=NA,minShape=NA,failureRegion=NA,threshold="otsu",numWindows=2){
	if(filename != ""){
		imageData=try(segmentImage(filename=filename,maxShape=maxShape,minShape=minShape,failureRegion=failureRegion,threshold=threshold,numWindows=numWindows))
	}else{
		imageData=try(segmentImage(image=image,maxShape=maxShape,minShape=minShape,failureRegion=failureRegion,threshold=threshold,numWindows=numWindows))
	}
	featuresObjects=imageData[[3]]
	image=imageData[[1]]
	segmentedImage=imageData[[2]]
	res=paintObjects(segmentedImage,image,col="green")
	nuctext=featuresObjects[,"index"]
	xy=featuresObjects[,c("m.cx","m.cy")]
	font = drawfont(weight=600, size=10)
	text=drawtext(res,xy=as.matrix(xy),labels=as.character(nuctext),font=font,col="yellow")
	l=list(text,featuresObjects)
}
