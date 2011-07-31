imageCompression <-
function(img,k){
	s1=svd(img[,,1],k,k)
	s2=svd(img[,,2],k,k)
	s3=svd(img[,,3],k,k)
	
	d1 <- diag(s1$d)
	d2 <- diag(s2$d)
	d3 <- diag(s3$d)
	u1=s1$u
	u2=s2$u
	u3=s3$u
	
	v1=s1$v
	v2=s2$v
	v3=s3$v
	
	C1=array(0,dim(img[,,1]))
	C2=array(0,dim(img[,,1])) 
	C3=array(0,dim(img[,,1]))
	for(j in 1:k){
		j
		C1=C1+d1[j,j] * u1[,j]%*% t(v1[,j])
	}
	for (j in 1:k){
		C2=C2+d2[j,j] * u2[,j] %*% t(v2[,j])
	}
	for (j in 1:k){
		C3=C3+d3[j,j] * u3[,j] %*% t(v3[,j])
	}
	img[,,1]=C1
	img[,,2]=C2
	img[,,3]=C3
	img
}

