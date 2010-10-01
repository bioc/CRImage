\name{processAperio}
\alias{processAperio}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cellularity Calculation of Aperio TX Scanner}
\description{
  Procession of Aperio TX Slides.
}
\usage{
processAperio(classifier = classifier, inputFolder = inputFolder, outputFolder = outputFolder, identifier = identifier, numSlides = numSlides, cancerIdentifier = cancerIdentifier,maxShape=NA,minShape=NA,failureRegion=NA,slideToProcess=NA,KS=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{classifier}{ The classifier.}
  \item{inputFolder}{ The path to the image folder. }
  \item{outputFolder}{ The path to the output folder. }
  \item{identifier}{ The identifier of the files ("Ss" or "Da") }
  \item{numSlides}{ The number of sections in the image. }
  \item{cancerIdentifier}{ The identifier of the cancer class}
\item{maxShape}{Maximum size of cell nuclei}
  \item{minShape}{ Minimum size of cell nuclei }
  \item{failureRegion}{ minimum size of failure regions}
  \item{slideToProcess}{Set this parameter if only a certain slide should be processed}
\item{KS}{Apply Kernel Smoother?}
}
\details{
  The function processes images of Aperio TX scanners. The images have to be saved in the CWS format. 
}
\value{
 Four folders are created in the output folder.
  \item{Files}{Cellularity values and cell numbers are saved in the file}
  \item{classifiedImage}{Subimages with labeled tumour and non tumour cells}
  \item{tumourDensity}{Cancer heatmaps for every subimage}
  \item{cellCoordinates}{Coordinates and cell class for every cell in the subimage}
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
\keyword{misc}
