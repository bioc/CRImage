\name{numberOfNeighbors}
\alias{numberOfNeighbors}
\title{
Internal function to calculate the number of neighbors of a cell.
}
\description{
Calculates the number of neighbors of a cell nuclei in a distance of 50 pixel.}
\usage{
numberOfNeighbors(img, cellCoordinates, allFeatures)
}
\arguments{
  \item{img}{
The image.}
  \item{cellCoordinates}{
Coordinates of the cells.}
  \item{allFeatures}{
Features of the cells.}
}
\details{
The function is an internal function.}
\value{
The number of neighbors for every cell.
}

\author{
Henrik Failmezger, failmezger@cip.ifi.lmu.de}
\examples{
#segment image
f = system.file('extdata' ,'exImg.jpg',package='CRImage')
segmentationValues=segmentImage(f,maxShape=800,minShape=40,failureRegion=2000)
image=segmentationValues[[1]]
segmentedImage=segmentationValues[[2]]
imageFeatures=segmentationValues[[3]]
}
\keyword{ misc }
