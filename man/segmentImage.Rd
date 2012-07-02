\name{segmentImage}
\alias{segmentImage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Segmentation of an image }
\description{
 The function segments cells or cell nuclei in the image.
}
\usage{
segmentImage(filename="",image=NA,maxShape=NA,minShape=NA,failureRegion=NA,threshold="otsu",window=2, colorCorrection=FALSE, classifyStructures=FALSE,pixelClassifier=NA,greyscaleImage=0,penClassifier=NULL,referenceHist=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{ A path to an image}
  \item{image}{ An 'image' object, if no filename is specified.}
  \item{maxShape}{Maximum size of cell nuclei}
  \item{minShape}{ Minimum size of cell nuclei }
  \item{failureRegion}{ minimum size of failure regions}
  \item{threshold}{Thresholding method, "otsu" or "phansalkar"}
\item{window}{Number of windows to use for thresholding.}
\item{colorCorrection}{deprecated}
\item{classifyStructures}{Segment structures in the image, if yes a pixel classifier has to be defined}
\item{pixelClassifier}{A SVM which classifies RGB color values in foreground and background.}
\item{greyscaleImage}{Channel of the RGB image, to use for thresholding, if 0 use a joined greyscale image.}
\item{penClassifier}{Classifier to exclude low quality images(will be part of next release)}
\item{referenceHist}{Color histogram of a reference image, that can be used to estimate the quality of the recent image (will be part of next release)}
}
\details{
The image is converted to greyscale and thresholded.
 Clutter is deleted using morphological operations. Clustered objects are separated using watershed algorithm.
 Segmented Cell nuclei, which exceed the maximum size are thresholded and segmented again. Cell nuclei which fall below the minimum size are
 deleted. Dark regions which exceed the parameter failureRegion are considered as artefacts and deleted.
 If the parameters are not defined, the operations will not be executed.
 Features are generated for every segmented object.
}
\value{
  A list is returned containing
  \item{image}{The original image}
  \item{segmented image}{The segmented image}
}
\references{EBImage, 'http://www.bioconductor.org/packages/release/bioc/html/EBImage.html' }
\author{ Henrik Failmezger, failmezger@cip.ifi.lmu.de }
\examples{
#segment image
#f = system.file('extdata' ,'exImg.jpg',package='CRImage')
#segmentationValues=segmentImage(f,maxShape=800,minShape=40,failureRegion=2000,threshold="otsu",numWindows=4)
#image=segmentationValues[[1]]
#segmentedImage=segmentationValues[[2]]
#imageFeatures=segmentationValues[[3]]
}
\keyword{misc}
