\name{calculateCellularity}
\alias{calculateCellularity}
\title{Calculation of tumour cellularity}
\description{
The function calculates the tumour cellularity of an image by counting tumour and non tumour cells.
}
\usage{
calculateCellularity(filename = "", image = NA, classifier, cancerIdentifier, KS=FALSE,maxShape=NA,minShape=NA,failureRegion=NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{A path to an image file.}
  \item{image}{If filename is undefined, an Image object }
  \item{classifier}{A SVM object, created with createClassifier or directly with the package e1071}
  \item{cancerIdentifier}{ A string which describes, how the cancer class is named.}
  \item{KS}{Apply kernel smoother?}
 \item{maxShape}{Maximum size of cell nuclei}
  \item{minShape}{ Minimum size of cell nuclei }
  \item{failureRegion}{ minimum size of failure regions}
}
\details{
  The method calculates tumour cellularity of an image. The cells of the image are classified and the cellularity is: numTumourCells/numPixel.
  Furthermore the number of cells of the different classes are counted. A heatmap of cellularity is created. The image is divided in 16 subwindows and cellularity is  calculated for every subwindow. Green in the heatmaps indicates strong cellularity, white low cellularity.
}
\value{
  A list containing
  \item{cellularity values }{a vector, the n first values indicate the n numbers of cells in the n classes, the n + 1th value indicates the tumour cellularity, The n + 2th value is the ratio of tumour cells by all cells}
  \item{cancerHeatmap }{Heatmap of cancer density}
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