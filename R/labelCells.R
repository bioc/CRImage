labelCells <-
function(img, segmentedImage,classes,classColours,nblocks=3,labeledPoints=NULL,filename=NULL,filenameImage=NULL,transformCoordinates=FALSE){
#Start positions for the hull
	oldX=NULL
	oldY=NULL
#
	options(stringsAsFactors = FALSE)
	f=computeFeatures.moment(segmentedImage)
	print(f)
	activateHull=FALSE
	convHull=c()
	xyCell=f[,c("m.cx","m.cy")]
	blockPosition=data.frame(stringsAsFactors=FALSE)
	xs=floor(dim(img)[1]/nblocks)
	ys=floor(dim(img)[2]/nblocks)
	blockCounter=1
	xPos=1
	yPos=1
	for (i in 1:nblocks){
		xPos=1
		for (j in 1:nblocks){
			blockPosition=rbind(blockPosition,c(blockCounter,xPos,(xPos+xs-1),yPos,(yPos+ys-1)))
			blockCounter=blockCounter+1
			xPos=xPos+(xs-1)
		}
		yPos=yPos+(ys-1)
	}
	colnames(blockPosition)=c("block","xs","xe","ys","ye")
	classCounter=1
	actClass=classes[classCounter]
	remove=FALSE
	action="add"
	description=c("Key: A=add point, D=delete point, H=Hull, C=switch class, Q=Exit, R=refresh, W=Write image to file ")
	actValues=paste("Class:",actClass,"Action:",action)
	actBlock=1;
	
#################
	if(is.null(labeledPoints)){
		print("No points labled.")
		labeledPoints= data.frame(t(rep(0,7)),stringsAsFactors=FALSE)
#you can change number_of_columns according to your need
		colnames(labeledPoints)=c("index","x","y","classCell","xLocal","yLocal","block")
	}else{
		labeledPoints= data.frame(as.matrix(labeledPoints),stringsAsFactors=FALSE) #transform to right coordinates
		if(transformCoordinates==TRUE){
			labeledPointsN=labeledPoints
			for(i in 1:nblocks){
				actPosition=blockPosition[blockPosition$block==i,]
				labeledPointsToTransformBlock=labeledPoints[labeledPoints$block==i,]
				xCoordLocal=abs(actPosition$ye-actPosition$ys)-as.numeric(labeledPointsToTransformBlock$xLocal)
				yCoordLocal=abs(actPosition$xe-actPosition$xs)-as.numeric(labeledPointsToTransformBlock$yLocal)
				labeledPoints[labeledPoints$block==i,"xLocal"]=yCoordLocal
				labeledPoints[labeledPoints$block==i,"yLocal"]=xCoordLocal
			}
		}
	}
#########
	labeledPointsToPaint=labeledPoints #these points are segmented nuclei and will painted
	labeledPointsNoNuclei=data.frame(stringsAsFactors=FALSE) #these nuclei will be marked
	if(dim(labeledPoints)[1]>1){
		if(length(which(as.numeric(labeledPoints[,1])==0))>0){
			labeledPointsToPaint=labeledPoints[-which(as.numeric(labeledPoints[,1])==0),]
			labeledPointsNoNuclei=labeledPoints[which(as.numeric(labeledPoints[,1])==0),]
		}
	}
	if(dim(labeledPointsToPaint)[1]>0){
		paintedNuclei=paintCells(segmentedImage,img,as.character(labeledPointsToPaint[,4]),as.numeric(labeledPointsToPaint[,1]),classes,colors=classColours)
	}
	actPosition=blockPosition[blockPosition$block==actBlock,]
	paintedNucleiN=paintedNuclei[actPosition$xs:actPosition$xe,actPosition$ys:actPosition$ye,]#y coordinates are traversed
	paintedNucleiN=aperm(paintedNucleiN,c(2,1,3))#transpose image to get correct view
#imgMatrix=imagematrix(paintedNucleiN)
#plot.imagematrix(imgMatrix)
	plotImage(paintedNucleiN)
	if(dim(labeledPointsNoNuclei)[1]>0){
		labeledPointsNoNucleiBlock=labeledPointsNoNuclei[labeledPointsNoNuclei$block==actBlock,]
		xCoord=as.numeric(labeledPointsNoNucleiBlock$xLocal)
		yCoord=as.numeric(labeledPointsNoNucleiBlock$yLocal)
		points(xCoord,yCoord,col="red",pch="x")
	}
	title(main=description,sub = actValues, col.sub="black",cex.main= 0.7)
################
	
	refresh=function(){
##refresh all nuclei
		trainingNuclei=segmentedImage
		trainingNuclei[,]=0;
		allIndices=c()
		
		labeledPointsToPaint=labeledPoints #these points are segmented nuclei and will painted
		labeledPointsNoNuclei=data.frame(stringsAsFactors=FALSE) #these nuclei will be marked
		if(dim(labeledPoints)[1]>1){
			if(length(which(as.numeric(labeledPoints[,1])==0))>0){
				labeledPointsToPaint=labeledPoints[-which(as.numeric(labeledPoints[,1])==0),]
				labeledPointsNoNuclei=labeledPoints[which(as.numeric(labeledPoints[,1])==0),]
			}
		}
		actPosition=blockPosition[blockPosition$block==actBlock,]
		print("Label nuclei")
		paintedNucleiSmall=paintCells(segmentedImage[actPosition$xs:actPosition$xe,actPosition$ys:actPosition$ye],img[actPosition$xs:actPosition$xe,actPosition$ys:actPosition$ye,],as.character(labeledPointsToPaint[,4]),as.numeric(labeledPointsToPaint[,1]),classes,colors=classColours)
		print("Nuclei labeled")
		paintedNuclei[actPosition$xs:actPosition$xe,actPosition$ys:actPosition$ye,]=paintedNucleiSmall
		paintedNuclei<<-paintedNuclei
		paintedNucleiSmall=aperm(paintedNucleiSmall,c(2,1,3))#transpose image to get correct view
#imgMatrix=imagematrix(paintedNucleiSmall)
#plot.imagematrix(imgMatrix)
		plotImage(paintedNucleiSmall)
#find block specific points and use local coordinates
		labeledPointsNoNucleiBlock=labeledPointsNoNuclei[labeledPointsNoNuclei$block==actBlock,]
		points(labeledPointsNoNucleiBlock$xLocal,labeledPointsNoNucleiBlock$yLocal,col="red",pch="x")
		title(main=description,sub = actValues, col.sub="black",cex.main= 0.7)
#save the file
		if(!is.null(filename)){
			write.table(labeledPoints,file=filename,sep="\t",row.names=FALSE,quote=FALSE)
		}
	}
	
	keydown <- function(key) {
		if (key == "r"){
			refresh()
		}
		if(key=="q"){
			if(!is.null(filename)){
				write.table(labeledPoints,file=filename,sep="\t",row.names=FALSE,quote=FALSE)
				print("Save")
			}
			print("Quit")
			return (invisible(1))
		}
		if(key=="w"){
			print("Write labeled image to file.")
			if(!is.null(filenameImage)){#label red crosses
				for (n in 1:dim(labeledPoints)[1]) {
					if(as.numeric(labeledPoints$index[n])==0){
						actNuclei=labeledPoints[n,]
						sizeNCrossU=3
						sizeNCrossD=3
						sizeNCrossR=3
						sizeNCrossL=3
						crossL=(as.numeric(actNuclei$x)-sizeNCrossL)
						crossR=(as.numeric(actNuclei$x)+sizeNCrossR)
						crossU=(as.numeric(actNuclei$y)-sizeNCrossU)
						crossD=(as.numeric(actNuclei$y)+sizeNCrossD)
						if((as.numeric(actNuclei$x)-sizeNCrossL)<1){
							crossL=1
						}
						if((as.numeric(actNuclei$y)-sizeNCrossU)<1){
							crossU=1
						}
						if((as.numeric(actNuclei$x)+sizeNCrossR)>dim(paintedNuclei)[1]){
							crossR=dim(paintedNuclei)[1]
						}
						if((as.numeric(actNuclei$y)+sizeNCrossD)>dim(paintedNuclei)[2]){
							crossD=dim(paintedNuclei)[2]
						}
						paintedNuclei[crossL:crossR,as.numeric(actNuclei$y),1]=1
						paintedNuclei[crossL:crossR,as.numeric(actNuclei$y),2]=0
						paintedNuclei[crossL:crossR,as.numeric(actNuclei$y),3]=0
						paintedNuclei[as.numeric(actNuclei$x),crossU:crossD,1]=1
						paintedNuclei[as.numeric(actNuclei$x),crossU:crossD,2]=0
						paintedNuclei[as.numeric(actNuclei$x),crossU:crossD,3]=0
					}
				}
				writeImage(paintedNuclei,filenameImage)
				message("Image saved")
			}else{
				message("You have to specify an image filename beforehand in order to save images.")
			}
		}
		if(key=="c"){
			classCounter<<-classCounter+1
			if(classCounter>length(classes)){classCounter<<-1}
#set old title to white
			title(sub = actValues, col.sub="white")
			actClass<<-classes[classCounter]
			actValues<<-paste("Class:",actClass,"Action:",action)
			print(paste("New class:",actClass))
#title(sub="test", col.sub="white") 
			title(sub = actValues, col.sub="black")
			NULL
		}
		if(key=="d"){
			print("Delete points")
			remove<<-TRUE
			title(sub = actValues, col.sub="white")
			actClass<<-classes[classCounter]
			action<<-"delete"
			actValues<<-paste("Class:",actClass,"Action:",action)
			title(sub = actValues, col.sub="black")
			NULL
		}
		if(key=="h"){
			
			if(activateHull==TRUE){
				message("Hull deactivated")
				title(sub = actValues, col.sub="white")
				action<<-""
				actValues<<-paste("Class:",actClass,"Action:",action)
				activateHull<<-FALSE
				title(sub = actValues, col.sub="black")
			}else{
				message("Hull activated")
				title(sub = actValues, col.sub="white")
				action<<-"Hull"
				actValues<<-paste("Class:",actClass,"Action:",action)
				activateHull<<-TRUE
				title(sub = actValues, col.sub="black")
			}
			NULL
		}
		if(key=="a"){
			print("Add points")
			remove<<-FALSE
			title(sub = actValues, col.sub="white")
			actClass<<-classes[classCounter]
			action<<-"add"
			actValues<<-paste("Class:",actClass,"Action:",action)
			title(sub = actValues, col.sub="black")
			NULL
		}
		if(key=="Right"){
			print("block")
			actBlock<<-actBlock+1
			if(actBlock>(nblocks*nblocks)){
				actBlock<<-1
			}
			refresh()
#actPosition=blockPosition[blockPosition$block==actBlock,]
#paintedNucleiN=paintedNuclei[actPosition$xs:actPosition$xe,actPosition$ys:actPosition$ye,]
#imgMatrix=imagematrix(paintedNucleiN)
#plot.imagematrix(imgMatrix)
#title(main=description,sub = actValues, col.sub="black",cex.main= 0.8)
			print(paste("Block number:",actBlock))
			NULL
		}
		if(key=="Left"){
			actBlock<<-actBlock-1
			if(actBlock<1){
				actBlock<<-(nblocks*nblocks)
			}
			refresh()
#actPosition=blockPosition[blockPosition$block==actBlock,]
#paintedNucleiN=paintedNuclei[actPosition$xs:actPosition$xe,actPosition$ys:actPosition$ye,]
#imgMatrix=imagematrix(paintedNucleiN)
#plot.imagematrix(imgMatrix)
#title(main=description,sub = actValues, col.sub="black",cex.main= 0.8)
			print(paste("Block number:",actBlock))
			NULL
		}
	}
	
	
#labeledPoints=data.frame()
#colnames(labeledPoints)=c("index","x","y","classCell")
#labeledPoints= data.frame(t(rep(0,4)))
#you can change number_of_columns according to your need
#	colnames(labeledPoints)=c("index","x","y","classCell")
#Assign column names (i have assumed number_of_columns as 6)
	
	
	
	
	
	dragmousedown <- function(buttons, x, y) {
		startx <- x
		starty <- y
		clicked.x <- try(grconvertX(x, from = "ndc",to = "user"))
		clicked.y <- try(grconvertY(y, from = "ndc",to = "user"))
		oldX<<-clicked.x
		oldY<<-clicked.y
		point=round(c(clicked.x,clicked.y))
		actPosition=blockPosition[blockPosition$block==actBlock,]
#transform local to global coordinates,coordinates are switched for image plots
		point[1]=actPosition$xs+point[1]
		point[2]=abs((actPosition$ye-actPosition$ys))-point[2]
		point[2]=actPosition$ys+point[2]
#
		if(point[1]>=1 & point[1]<dim(segmentedImage)[1]&point[2]>=1 & point[2]<dim(segmentedImage)[2]){
			index=segmentedImage[point[1],point[2]]
			if(!activateHull){
				if(remove==FALSE){
					print("add point")
					print(paste("X:",round(point[1]),"Y:",round(point[2]),"Class:",actClass,"Nucleus:",index))
#check if index already exists
#nucleiIndexToAdd=which(labeledPoints$index==index)
#if exists delete it
#if(length(nucleiIndexToAdd)>0 & index != 0){
#	print("Replace label")
#	labeledPoints<<-labeledPoints[-nucleiIndexToAdd,]
#}
					labeledPoints<<-rbind(labeledPoints,c(index,point[1],point[2],actClass,round(clicked.x),round(clicked.y),actBlock))
				}else{
					print("remove point")
					print(paste("X:",round(point[1]),"Y:",round(point[2]),"Class:",actClass,"Nucleus:",index))
					if(index>0){
						nucleiIndexToRemove=which(labeledPoints$index==index)
#check if this index exists
						if(length(nucleiIndexToRemove)>0){
							labeledPoints<<-labeledPoints[-nucleiIndexToRemove,]
						}
#if data frame empty do new initialisation
						if(dim(labeledPoints)[1]==0){
							labeledPoints= data.frame(t(rep(0,7)),stringsAsFactors=FALSE)
#you can change number_of_columns according to your need
							colnames(labeledPoints)=c("index","x","y","classCell","xLocal","yLocal","block")
							labeledPoints<<-labeledPoints
						}
					}
				}
			}
		}else{
			message("Point outside device")
		}
		NULL
#devset()
#	usr <<- par("usr")	
	}
#########################
	
	dragmousemove=function(buttons, x, y) {
		if(activateHull){
			clicked.x <- grconvertX(x, from = "ndc",to = "user")
			clicked.y <- grconvertY(y, from = "ndc",to = "user")
			
#print(paste("x:",clicked.x))
#print(paste("y:",clicked.y))
#print(paste("oldX:",oldX))
#print(paste("oldY:",oldY))
#print(paste("dist",sqrt((clicked.x-oldX)^2+(clicked.y-oldY)^2)))
			if(!is.null(oldX) & !is.null(oldY)){
				if(sqrt((clicked.x-oldX)^2+(clicked.y-oldY)^2)>0.5){
					"draw"
#print(paste("dist",sqrt((clicked.x-oldX)^2+(clicked.y-oldY)^2)))
					
					lines(c(oldX,clicked.x),c(oldY,clicked.y),col=classColours[which(classes==actClass)])
					convHull<<-rbind(convHull,c(clicked.x,clicked.y))
					oldX<<-clicked.x
					oldY<<-clicked.y
				}
			}
		}
		NULL
	}
	dragmouseup <- function(buttons, x, y) {
		if(activateHull){
			if(!is.null(dim(convHull)[1])){
				actPosition=blockPosition[blockPosition$block==actBlock,]
				
				convHull[,2]=abs((actPosition$ye-actPosition$ys)-convHull[,2])
				convHull[,2]=actPosition$ys+convHull[,2]
				convHull[,1]=actPosition$xs+convHull[,1]
				
				convHullS=convHull[sample(1:dim(convHull)[1],dim(convHull)[1],replace=F),]
#
				
				actPointsXYBlockIndex=which(xyCell[,1]<actPosition$xe & xyCell[,1]>=actPosition$xs & xyCell[,2]<actPosition$ye & xyCell[,2]>=actPosition$ys)
				actPointsXYBlock=xyCell[actPointsXYBlockIndex,]
#transform local to global coordinates,coordinates are switched for image plots
#
#tr<-tri.mesh(convHullS[,1],convHullS[,2],"strip")
				tr=chull(convHullS[,1],convHullS[,2])
#pointsInHull=in.convex.hull(tr,actPointsXYBlock[,1],actPointsXYBlock[,2])
				pointsInHull=in.chull(actPointsXYBlock[,1],actPointsXYBlock[,2],convHullS[tr,1],convHullS[tr,2])
				indexPointsInHull=actPointsXYBlockIndex[which(pointsInHull==1)]
				
#print(pointsInHull)
#print(index)
#print(replicate(length(index),actClass))
				newPoints=cbind(indexPointsInHull,actPointsXYBlock[pointsInHull,1],actPointsXYBlock[pointsInHull,2],replicate(length(indexPointsInHull),actClass),actPointsXYBlock[pointsInHull,1],actPointsXYBlock[pointsInHull,2],replicate(length(indexPointsInHull),actBlock))
#######
#check if index already exists
				existingIndices=intersect(labeledPoints$index,indexPointsInHull)
#if exists delete it
				if(length(existingIndices)>0){
					message("Replace label")
					labeledPoints<<-labeledPoints[-existingIndices,]
				}
				
######
				colnames(newPoints)=c("index","x","y","classCell","xLocal","yLocal","block")
#print(newPoints)
#print(labeledPoints)
				labeledPoints<<-rbind(labeledPoints,newPoints)
#print(labeledPoints)
				refresh()
				clicked.x <- grconvertX(x, from = "ndc",to = "user")
				clicked.y <- grconvertY(y, from = "ndc",to = "user")
				oldX<<-NULL
				oldY<<-NULL
				convHull<<-c()
			}
		}
#usr <<- par("usr")
		NULL
	}
	
	setGraphicsEventHandlers(prompt="Click to label cells, hit q to quit",onMouseMove = dragmousemove,onMouseUp = dragmouseup,onMouseDown = dragmousedown,onKeybd = keydown)	
#eventEnv <- getGraphicsEventEnv()
	getGraphicsEvent()
#labeledPoints=labeledPoints[-which(labeledPoints[,1]==0),]
	return(labeledPoints)
}

