hist3d <-
function(imgS,bins,whitePixelMask){
    HN=imgS
	RChannel=HN[,,1]
#RChannel=round(RChannel*bins*(1/(256-1)))
	GChannel=HN[,,2]
#	GChannel=round(GChannel*bins*(1/(256-1)))
	BChannel=HN[,,3]
#	BChannel=round(BChannel*bins*(1/(256-1)))
	t_count=floor(RChannel[whitePixelMask == FALSE]*bins)+256*(floor(bins*BChannel[whitePixelMask == FALSE]))+256*256*(floor(bins*GChannel[whitePixelMask == FALSE]))
	colorCounts=table(as.vector(t_count))
	histogram=array(0,c(bins+1,bins+1,bins+1))
	for(i in 1:length(colorCounts)){
		value1=as.numeric(names(colorCounts)[i])%%256
		value2=floor((as.numeric(names(colorCounts)[i])%%(256*256))/256)
		value3=floor((as.numeric(names(colorCounts)[i])/(256*256)))
		if((value1) >(bins)){value1=bins}
		if((value2) >(bins)){value2=bins}
		if((value3) >(bins)){value3=bins}
		histogram[value1+1,value2+1,value3+1]=colorCounts[i]
	}
	histogram=rbind(as.vector(histogram[,,1]),as.vector(histogram[,,2]),as.vector(histogram[,,3]))
}
