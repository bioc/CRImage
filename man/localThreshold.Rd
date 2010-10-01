\name{localThreshold}
\alias{localThreshold}
\title{
Internal function to do local thresholding.
}
\description{
The function calculates a local threshold of an image using calculateThreshold.}
\usage{
localThreshold(imgG, img)
}
\arguments{
  \item{imgG}{
Greyscale image.}
  \item{img}{
Image with colors.}
}
\details{
The function is an internal function.
}
\value{
The thresholded image.
}

\author{
Henrik Failmezger, <email: failmezger@cip.ifi.lmu.de>}
\examples{
#segment image
f = system.file('extdata' ,'exImg.jpg',package='CRImage')
segmentationValues=segmentImage(f,maxShape=800,minShape=40,failureRegion=2000)
image=segmentationValues[[1]]
segmentedImage=segmentationValues[[2]]
imageFeatures=segmentationValues[[3]]
}
\keyword{ misc }
