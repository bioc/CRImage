convertHSVToRGB <-
function(imgHSV){
	
	h=as.vector(imgHSV[,,1])
	s=as.vector(imgHSV[,,2])
	v=as.vector(imgHSV[,,3])
	hi = floor(h/60);
	
	f=(h/60)-hi
	
	p=v*(1-s)
	q=v*(1-s*f)
	t=v*(1-s*(1-f))
	
	
	r=rep(0,length(h))
	g=rep(0,length(h))
	b=rep(0,length(h))
	
	ind=which(hi==0)
	r[ind]=v[ind]
	g[ind]=t[ind]
	b[ind]=p[ind]
	
	ind=which(hi==6)
	r[ind]=v[ind]
	g[ind]=t[ind]
	b[ind]=p[ind]
	
	ind=which(hi==1)
	r[ind]=q[ind]
	g[ind]=v[ind]
	b[ind]=p[ind]
	
	ind=which(hi==2)
	r[ind]=p[ind]
	g[ind]=v[ind]
	b[ind]=t[ind]
	
	
	ind=which(hi==3)
	r[ind]=p[ind]
	g[ind]=q[ind]
	b[ind]=v[ind]
	
	ind=which(hi==4)
	r[ind]=t[ind]
	g[ind]=p[ind]
	b[ind]=v[ind]
	
	
	ind=which(hi==5)
	r[ind]=v[ind]
	g[ind]=p[ind]
	b[ind]=q[ind]
	
	imgHSV[,,1]=array(r,dim(imgHSV[,,1]))
	imgHSV[,,2]=array(g,dim(imgHSV[,,1]))
	imgHSV[,,3]=array(b,dim(imgHSV[,,1]))
	imgHSV
}

