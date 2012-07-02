convertLABToRGBOld <-
function(LAB){
	labToLMS1=matrix(c(1,1,1,1,1,-2,1,-1,0),c(3,3))
	labToLMS2=matrix(c(sqrt(3)/3,0,0,0,sqrt(6)/6,0,0,0,sqrt(2)/2),c(3,3))
	LMS=labToLMS1 %*% labToLMS2 %*% LAB
	LMS=LMS
	lmsToRGB=matrix(c(4.4679,-1.2186,0.0497,-3.5873,2.3809,-0.2439,0.1193,-0.1624,1.2045),c(3,3))
	RGB=lmsToRGB %*% LMS
	RGB
}
