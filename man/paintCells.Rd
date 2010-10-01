\name{paintCells}
\alias{paintCells}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Coloring of classified cells}
\description{
internal function
}
\usage{
paintCells(imgWT, img, classes, index, classValues)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{imgWT}{ segmented image }
  \item{img}{ the image }
  \item{classes}{ the classes }
  \item{index}{ index }
  \item{classValues}{ the class values }
}
\details{
internal function}

\author{ Henrik Failmezger,failmezger@cip.ifi.lmu.de }
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
\keyword{misc}
