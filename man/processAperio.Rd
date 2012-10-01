\name{processAperio}
\alias{processAperio}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cellularity Calculation of Aperio TX Scanner}
\description{
  Procession of Aperio TX Slides.
}
\usage{
processAperio(classifier=classifier,inputFolder=inputFolder,outputFolder=outputFolder,identifier=identifier,numSlides=numSlides,cancerIdentifier=cancerIdentifier,classOther = NA,maxShape=800,minShape=40,failureRegion=2000,slideToProcess=NA,KS=TRUE,colors=c(),classesToExclude=c(),threshold="otsu",numWindows=2,colorCorrection=FALSE,classifyStructures=FALSE,ksToExclude=c(),pixelClassifier=NA,densityToExclude=c(),numDensityWindows=32,resizeFactor=4,plotCellTypeDensity=TRUE,greyscaleImage=0,penClassifier=NULL,referenceHist=NULL,fontSize = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{classifier}{ The classifier.}
  \item{inputFolder}{ The path to the image folder. }
  \item{outputFolder}{ The path to the output folder. }
  \item{identifier}{ The identifier of the files ("Ss" or "Da") }
  \item{numSlides}{ The number of sections in the image. }
  \item{cancerIdentifier}{ The identifier of the cancer class}
\item{classOther}{deprecated}
\item{maxShape}{Maximum size of cell nuclei}
  \item{minShape}{ Minimum size of cell nuclei }
  \item{failureRegion}{ minimum size of failure regions}
  \item{slideToProcess}{Set this parameter if only a certain slide should be processed}
\item{KS}{Apply Kernel Smoother?}
 \item{colors}{Colors to paint the classes}
\item{classesToExclude}{Which class should be excluded?}
\item{threshold}{Which thresholding method should be used, "otsu" or "phansalkar" possible}
\item{numWindows}{Number of windows to use for thresholding.}
\item{colorCorrection}{deprecated}
\item{classifyStructures}{Use hierarchical classification. If yes a pixel classifier has to be defined. }
\item{ksToExclude}{These classes are excluded from kernel smoothing.}
\item{pixelClassifier}{A SVM to classify pixel based on their color values. Needed if hierarchical classification should be applied.}
\item{densityToExclude}{This class is excluded from cellularity calculation.}
\item{numDensityWindows}{Number of windows for the density plot.}
\item{resizeFactor}{Specifies the size of the cell density image. If this variable is not defined, the size of the thumbnail is used for the cell density image, else the size is calculated by size(thumbnail)*resizeFactor. The thumbnail is the small overview image, created by the Aperio software.}
\item{plotCellTypeDensity}{Plot the density of different cell types?}
\item{greyscaleImage}{Color channel of the RGB image that should be used for thresholding}
\item{penClassifier}{Classifier to exclude low quality images (will be part of next release)}
\item{referenceHist}{Colour Histogram of a reference image that can be used to calculate the quality of the recent image. (will be part of next release)}
\item{fontSize}{will be part of next release}
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
\item{resizeFactor}{Size of the cellularity density image, calculated by size(thumbnail) * resizeFactor. Whereas the thumbnail is the small overview image produced by Aperio.}
}
\author{ Henrik Failmezger, failmezger@mpipz.mpg.de}


\examples{
#t = system.file("extdata", "trainingData.txt", package="CRImage")
#read training data
#trainingData=read.table(t,header=TRUE)
#create classifier
#classifier=createClassifier(trainingData,topo=FALSE)[[1]]
#classify aperio
#f = system.file("extdata",  package="CRImage")
#f=file.path(f,"8905")
#dir.create("AperiOutput")
#takes long time!

f = system.file("extdata",  package="CRImage")
fc=file.path(f,"testClassifier")
load(fc)
fp=file.path(f,"pixelClassifier")
load(fp)
pixelClassifier=model
pathToImage=file.path(f,"8905")

pathToOutput="" #specify an output folder here

#processAperio(classifier=classifier,inputFolder=pathToImage,outputFolder=pathToOutput,identifier="Da",numSlides=1,cancerIdentifier="c",maxShape=800,minShape=40,failureRegion=2000,slideToProcess=1,KS=FALSE,colors=c("white","green","blue","brown"),classesToExclude=c('a','l','nl'),threshold="otsu",ksToExclude=c('l','nl'),pixelClassifier=pixelClassifier,classifyStructures=TRUE,densityToExclude=c('a'),numDensityWindows=8,resizeFactor=1.5,plotCellTypeDensity=FALSE,greyscaleImage=1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
