\name{segmentCytoplasma}
\alias{segmentCytoplasma}
\title{
Internal function to segment cell cytoplasma
}
\description{
The functions segments the cytoplasm of a cell.}
\usage{
segmentCytoplasma(img, imgW, indexWhitePixel, imgG, index, hF)
}
\arguments{
  \item{img}{
The image.}
  \item{imgW}{
The segmented image.}
  \item{indexWhitePixel}{
Index of white pixels.}
  \item{imgG}{
The greyscale image.}
  \item{index}{
THe index of the cells.}
  \item{hF}{
The hull features of the cells}
}
\details{
The function is an internal function.}
\value{
The segmented cytoplasma
}

\author{
Henrik Failmezger, <email: failmezger@cip.ifi.lmu.de>
}
\examples{
#segment image
f = system.file('extdata' ,'exImg.jpg',package='CRImage')
segmentationValues=segmentImage(f,maxShape=800,minShape=40,failureRegion=2000)
image=segmentationValues[[1]]
segmentedImage=segmentationValues[[2]]
imageFeatures=segmentationValues[[3]]
}
\keyword{ misc }
