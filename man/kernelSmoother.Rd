\name{kernelSmoother}
\alias{kernelSmoother}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ A kernel smoother to classify images }
\description{
  internal function
}
\usage{
kernelSmoother(predictedClasses, cellCoordinates, segmentedImage, cancerIdentifier, indexCells)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{predictedClasses}{ the predicted classes}
  \item{cellCoordinates}{ the coordinates of the cells}
  \item{segmentedImage}{ the segmented image}
  \item{cancerIdentifier}{ the label of the cancer identifier }
  \item{indexCells}{ the index of the cells }
}
\details{
internal function
}
\value{
the classes calculated by the kernel smoother
}
\author{ Henrik Failmezger, failmezger@cip.ifi.lmu.de }
\examples{
t = system.file("extdata", "trainingData.txt", package="CRImage")
#read training data
trainingData=read.table(t,header=TRUE)
#create classifier
classifier=createClassifier(trainingData,topo=FALSE)[[1]]
#classify cells
f = system.file("extdata", "exImg.jpg", package="CRImage")
classesValues=classifyCells(classifier,filename=f,KS=TRUE,,maxShape=800,minShape=40,failureRegion=2000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
