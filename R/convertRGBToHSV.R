convertRGBToHSV <-
function(img){
	imgHSV=imageData(img)
	imgR=as.vector(img[,,1])
	imgG=as.vector(img[,,2])
	imgB=as.vector(img[,,3])
	
	Max_rgb=apply(cbind(imgR,imgG,imgB),1,max)
	Min_rgb=apply(cbind(imgR,imgG,imgB),1,min)
	
	v = Max_rgb
	
	h=rep(0,length(v))
	s=v-Min_rgb
	
	
	
	k=which(imgR==v)
#h[k]=60*((imgR[k]-imgB[k])/s[k])
	h[k]=60*((imgG[k]-imgB[k])/s[k])
	
	k=which(imgG==v)
#	h[k]=60*(2+(imgG[k]-imgB[k])/s[k])
	h[k]=60*(2+(imgB[k]-imgR[k])/s[k])
	
	k=which(imgB==v)
	h[k]=60*(4+(imgR[k]-imgG[k])/s[k])
	
	k=which(Max_rgb==Min_rgb)
	h[k]=0
	
	
	k=which(h < 0)
	h[k]=h[k]+360
	
	k=which(Max_rgb==0)
	s[k]=0
	
	s=s/Max_rgb
	imgHSV[,,1]=h
	imgHSV[,,2]=s
	imgHSV[,,3]=v
	imgHSV
}
