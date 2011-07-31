\name{calculateMeanStdTarget}
\alias{calculateMeanStdTarget}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates Mean and Standard deviation of an image
}
\description{
Mean and SD calculation}
\usage{
calculateMeanStdTarget(imgT)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{imgT}{
the Image to calculate.
}
}
\details{
Mean and SD}
\value{
Vector with mean and standard deviation.
}

\author{
Henrik Failmezger, failmezger@cip.ifi.lmu.de
}


\examples{
	#read the target image
	f1= system.file("extdata", "exImg2.jpg", package="CRImage")
	targetImage=readImage(f1)
	#read the image whose color values should be adapted
	f2= system.file("extdata", "exImg3.jpg", package="CRImage")
	imgToConvert=readImage(f2)
	#calculate mean and standard deviation of target color channels
	mst=calculateMeanStdTarget(targetImage)
	# create a white pixel mask
	whitePixelMask=imgToConvert[,,1]>0.85 & imgToConvert[,,2]>0.85 & imgToConvert[,,3]>0.85
	#adapt color channels of image
imgCorrected=colorCorrection(imgToConvert,mst,whitePixelMask)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
\keyword{ misc }% __ONLY ONE__ keyword per line
