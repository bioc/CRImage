parseFinalScan <-
function(file){
	d.temp=scan(file,what="character",sep="\n")
	
	positionMatrix=data.frame(stringsAsFactors=FALSE)
	level=FALSE
	for(line in d.temp){
		if(substr(line,1,12)=="tDescription"){
			description=strsplit(line," ")[[1]]
			size=description[6]
			size=strsplit(description[6],"x")
			width=size[1]
			height=size[2]
			row=c("size",width[[1]],height[[1]])
			positionMatrix=row
		}
		if(substr(line,1,8)=="[Level0]"){
			level=TRUE
		}
		if(substr(line,1,1)=="[" && level==TRUE){
			row=c()
			splittedLine1=strsplit(line,"\\[")[[1]]
			splittedLine2=strsplit(splittedLine1[2],"\\]")[[1]]
			name=splittedLine2
		}
		if(substr(line,1,1)=="x" && level==TRUE){
			splittedLine1=strsplit(line,"=")[[1]]
			xValue=splittedLine1[[2]]
		}
		if(substr(line,1,1)=="y" && level==TRUE){
			splittedLine1=strsplit(line,"=")[[1]]
			yValue=splittedLine1[[2]]
			row=c(name,xValue,yValue)
			positionMatrix=rbind(positionMatrix,row)
		}
		
	}
	positionMatrix
}

