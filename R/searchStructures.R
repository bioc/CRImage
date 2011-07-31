searchStructures <-
function(img,imgW,imgSt,index,hF){
	
	smallImageX=250
	smallImageY=250
		
	minX=1
	minY=1


	xscale=dim(img)[1]-1
	#xscale=hF[,"g.x"]-minX
	#xscale=range(xscale)[2]-range(xscale)[1]

	
	#yscale=hF[,"g.y"]-minY
	#yscale=range(yscale)[2]-range(yscale)[1]
	yscale=dim(img)[1]-1

	featuresStructure=hullFeatures(imageData(imgSt))[,"g.s"]
	
	f=function(index,imgWt,imgBCsDt){
		
		x=round(as.vector(hF[index,"g.x"]))
		#calculate position in small Structure Image
		xS=x-1
		xS=(xS/(xscale))*(smallImageX-1)
		xS=xS+1
		y=round(as.vector(hF[index,"g.y"]))
		yS=y-1
		yS=(yS/yscale)*(smallImageY-1)
		yS=yS+1
		xS=round(xS)
		yS=round(yS)
		if(xS>=1 && yS>=1 && xS<=dim(imgSt)[1] &&  yS<=dim(imgSt)[2]){
			segSt=imgSt[xS,yS]
			if(segSt >0){
				#size of the structure is the feature
				segSt=featuresStructure[segSt]
			}else{
				0
			}
		}else{
			0
		}
	}
	imgBCsD=NULL
	structures=sapply(index,f, imgWt=imgW,imgBCsD=as.vector(imgBCsD),simplify=TRUE)
	
	structures

}

