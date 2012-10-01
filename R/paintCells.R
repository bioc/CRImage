paintCells <-
function(imgWT,img,classes,index,classValues,colors=c()){
	if(length(classValues)<=10){
		message("paint cells")
		indexClass=data.frame(index,classes,stringsAsFactors=FALSE)
		indexClass=rep("n",max(imgWT))
		indexClass[index]=classes
		indexClass=c(indexClass,"n")
		imgOld=img
		if(length(colors)==0){
			cols=c("white","green","red","yellow","blue","orange","brown","cyan","gray","purple","violet")
		}else{
			if(length(colors)<length(classValues)){
				message("Not enough colors specified for the class values. Standard coloring is used.")
				cols=c("white","green","red","yellow","blue","orange","brown","cyan","gray","purple","violet")
			}else{
				cols=colors
			}
		}
		counter=1
		for (c in classValues){
			imgTC=imgWT
			imgTC[imgTC==0]=length(indexClass)
			a=array(indexClass[imgTC] !=c ,dim(imgWT))
			imgTC[a]=0
			imgTC[imgTC==length(indexClass)]=0
			imgOld = paintObjects(imgTC,imgOld,col=cols[counter])
			counter=counter+1
		}
		
		imgOld
	}else{
		message("Too many classes to paint. Do not use more than 10 classes.")
		img
	}
}
