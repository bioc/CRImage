classifyPen <-
function(imgS,penClassifier,whitePixelMask){
	message("classify for pen containment")
	histImage=hist3d(imgS,10,whitePixelMask)
	histImage=histImage/sum(histImage)
	naIndex=which(is.na(histImage))
	if(length(naIndex)==0){
		label=predict(penClassifier,t(as.vector(histImage)^(1/4)))
	}else{
		label="n"
	}
	label
}
