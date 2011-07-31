Phansalkar_threshold <-
function(allGreyValues){
	m=mean(allGreyValues)
	s=sd(allGreyValues)
	q=10
	p=2
	k=0.25
	R=0.5
	t=m*(1+p*exp(-q*m)+k*((s/R)-1))
	t
}

