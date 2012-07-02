calculateOtsu <-
function(allGreyValues){
	message("Otsu threshold")
	gTrans=allGreyValues*255
	
	gTrans=round(gTrans)
	"h=hist(as.vector(gTrans))"
	t=tabulate(gTrans)
	t=t+0.001
	numPixel=length(allGreyValues)
	ni=t
	pi=ni/numPixel
	uT=function(p,L){
		res=0
		for (i in 1:L){
			res = res+i*p[i]
		}
		res
	}
	wK = function(p,k){
		res=0
		for (i in 1:k){
			res = res+p[i]
		}
		res
	}
	uK = function(p,k){
		res=0
		for (i in 1:k){
			res = res+i*p[i]
		}
		res
	}
	maxT=function(p,L,k){
		(uT(p,L)*wK(p,k)-uK(p,k))^2/(wK(p,k)*(1-wK(p,k)))
	}
	if(length(t)>1){
		o=optimize(maxT,c(1,length(t)),maximum=TRUE,L=length(t),p=pi)
		threshold=o$maximum/255
	}else{
		message("No values to optimize.")
		threshold=0
	}
	
}
