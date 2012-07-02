getImageDistance <-
function(imgS,referenceHist,whitePixelMask){
	message("calculate distance to ref")
	histImage=hist3d(imgS,10,whitePixelMask)
	histImage=histImage/sum(histImage)
	dist=sum((referenceHist-histImage)^2)
	dist
}
