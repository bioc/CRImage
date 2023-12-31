segmentCytoplasma <-
function(img,imgW,indexWhitePixel,imgG,index,hF){
	colorMode(img)="gray"
	imgGn=as.vector(img[,,1])
	imgGn[img[,,1]==2]=-1
	imgGn[imgGn==1]=-1
	imgGn[imgW>1]=-1
	imgGn=imgGn[-indexWhitePixel]
	t=calculateOtsu(as.vector(imgGn[imgGn !=-1]))
	imgBC=imgG
	imgBC[img[,,1]<=t]=-1
	imgBC[imgBC != -1]=0
	imgBC[imgBC == -1]=1
	imgBCo=opening(imgBC,makeBrush(5, shape='disc'))
	imgBCs=bwlabel(imgBCo)
	regions=unique(as.vector(imgBCs))
	if(length(regions)>0){
		imgBCsD=imageData(imgBCs)
		colorMode(img)="color"	
		"imgCyto=paintObjects(imgBCs,img,col=\"yellow\")"
		featuresCytoplasma=computeFeatures.shape(imageData(imgBCs))[,"s.area"]
		imgBCs=imageData(imgBCs)
		f=function(index,imgWt,imgBCsDt){
			x=round(as.vector(hF[index,"m.cx"]))
			y=round(as.vector(hF[index,"m.cy"]))
			if(x>0 & y>0){
				segCyto=imgBCs[x,y]
				if(segCyto >0){
					s=featuresCytoplasma[segCyto]
				}else{
					0
				}
			}else{
				0
			}
		}
	sizeCytoplasma=sapply(index,f, imgWt=imgW,imgBCsD=as.vector(imgBCsD),simplify=TRUE)
	}else{
	message("No cytoplasma found")
	sizeCytoplasma=rep(0,length(index))
	}
	sizeCytoplasma
}
