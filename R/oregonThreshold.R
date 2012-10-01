oregonThreshold <-
function(imgG){
	rPixelValues=as.vector(imgG)
	thrsh=mean(rPixelValues)/2
	removed=TRUE
	while(removed==TRUE){
		pixelSmaller=which(rPixelValues<thrsh)
		if(length(pixelSmaller)==0){
			removed==FALSE
			break
		}else{	
			rPixelValues=rPixelValues[-pixelSmaller]
			thrsh=mean(rPixelValues)/2
		}
	}
	thrsh
}
