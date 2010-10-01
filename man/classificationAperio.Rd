\name{classificationAperio}
\alias{classificationAperio}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Classification of ScanScopeTX Slides.}
\description{
The slides are segmented and classified. 
}
\usage{
classificationAperio(fileLocation, filename, pathToOutputFolderImgDir, classifier, pathToOutputFolderImgDirFiles, pathToOutputFolderImgDirImages, pathToOutputFolderImgDirCellDensity, blockSlice, sliceColors, sizeO, index, cancerIdentifier,maxShape=NA,minShape=NA,failureRegion=NA,KS=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fileLocation}{ location of the image file }
  \item{filename}{ name of the image file }
  \item{pathToOutputFolderImgDir}{ path where the classified images are saved }
  \item{classifier}{ the classifier object }
  \item{pathToOutputFolderImgDirFiles}{ the path where the files with the cellularity values should be saved.}
  \item{pathToOutputFolderImgDirImages}{ the path where the classified images are saved.}
  \item{pathToOutputFolderImgDirCellDensity}{ the path where the cancer heatmaps are saved }
  \item{blockSlice}{which slide does the image have }
  \item{sliceColors}{colors to label the classes}
  \item{sizeO}{ size of the image }
  \item{index}{ index of the image }
  \item{cancerIdentifier}{ label of the cancer identifier }
  \item{maxShape}{Maximum size of cell nuclei}
  \item{minShape}{ Minimum size of cell nuclei }
  \item{failureRegion}{ minimum size of failure regions}
  \item{KS}{Apply KernelSmoother?}
}
\details{
The function is an internal function which is used by processAperio. 
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
