SauvolaThreshold <-
function(allGreyValues){
	m=mean(allGreyValues)
	s=sd(allGreyValues)
	k=0.25
	R=0.5
	t=m*(1+k*((s/R)-1))
	t
}

