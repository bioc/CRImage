\name{findSlices}
\alias{findSlices}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{calculates the coordinates of the subimages }
\description{
internal function 
}
\usage{
findSlices(imgFolder, pathToOutputFolder, numSlides)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{imgFolder}{ folder to the subimages }
  \item{pathToOutputFolder}{ path to the output folder}
  \item{numSlides}{ number of sections in the image }
}
\details{
internal function
}
\value{
the sections for every image
}

\author{ Henrik Failmezger, failmezger@cip.ifi.lmu.de }

\examples{
t = system.file("extdata", "trainingData.txt", package="CRImage")
#read training data
trainingData=read.table(t,header=TRUE)
#create classifier
classifier=createClassifier(trainingData,topo=FALSE)[[1]]
#classify aperio
f = system.file("extdata",  package="CRImage")
f=file.path(f,"8905")
dir.create("AperiOutput")
#takes long time!
#processAperio(classifier=classifier,inputFolder=f,outputFolder="AperiOutput",identifier="Da",numSlides=2,cancerIdentifier="1",maxShape=800,minShape=40,failureRegion=2000)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
