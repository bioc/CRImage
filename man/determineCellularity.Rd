\name{determineCellularity}
\alias{determineCellularity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Determination of cellularity }
\description{
The function is an internal function which calculates cellularity.
}
\usage{
determineCellularity(classes, classifiedCells, dimImg, img, imgW, indexWhitePixel, cancerIdentifier,classValues)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{classes}{ the classes }
  \item{classifiedCells}{ the classified cells }
  \item{dimImg}{ dimension of the image }
  \item{img}{ the image }
  \item{imgW}{ the segmented image}
  \item{indexWhitePixel}{ the index of the white pixels }
  \item{cancerIdentifier}{ the label of the cancer identifier}
  \item{classValues}{ the class values }
}
\details{
The function is an internal function.
}
\value{
The calculated cellularity
}
\author{ Henrik Failmezger, failmezger@cip.ifi.lmu.de }
\examples{
t = system.file("extdata", "trainingData.txt", package="CRImage")
#read training data
trainingData=read.table(t,header=TRUE)
#create classifier
classifier=createClassifier(trainingData,topo=FALSE)[[1]]
#calculation of cellularity
f = system.file("extdata", "exImg.jpg", package="CRImage")
exImg=readImage(f)
cellularityValues=calculateCellularity(f,classifier=classifier,cancerIdentifier="1",maxShape=800,minShape=40,failureRegion=2000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
