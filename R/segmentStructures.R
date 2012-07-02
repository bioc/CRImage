segmentStructures <-
function(img,pixelClassifier){
	meanStdTarget=rbind(c(78.282154,300.371584) ,c(9.694320,-10.856946), c(2.081496,3.614328))
	message("Segment Structures")
	f = makeBrush(5, shape='disc', step=FALSE)
	f = f/sum(f)
	img=resize(img,250,250)
######TODELETE#####
	imgT=img
	imgCor=colorCorrectionOld(img,meanStdTarget)
	img=imgCor[[1]]
	imgC=img
	imgG=filter2(imgC, f)
	colorValues=cbind(as.vector(imgG[,,1]),as.vector(imgG[,,2]),as.vector(imgG[,,3]))
	predictedClasses=predict(pixelClassifier,colorValues)
	imgB=array(as.numeric(as.character(predictedClasses)),dim(imgC)[1:2])
	imgS=bwlabel(imgB)
	numSeq=tabulate(imageData(imgS)+1)
	imgSdN=imageData(imgS)+1
#set the failure region
	a=array(numSeq[imgSdN]<80,dim(imgSdN))
	imgSdN[a]=1
	imgSdN=imgSdN-1
	imgS=imgSdN
	imgP=paintObjects(imgS,imgT,col=c("lightgreen"))
#	allStructures=unique(imgS)
#	for (s in allStructures){
#meanS=mean(hueValue[which(imgS==s)])
#		print(meanS)
#if(meanS>=290){
#			imgS[imgS==s]=0
#		}else{
#			meanStructures=c(meanStructures,meanS)
	
#		}
#}
	l=list(imgS,imgB,imgG,imgP,imgC,imgT)
}
