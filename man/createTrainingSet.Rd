\name{createTrainingSet}
\alias{createTrainingSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Construction of a training set }
\description{
  Creates a training set for cell classification.
}
\usage{
createTrainingSet(filename = "", image = NA,maxShape = NA,minShape = NA,failureRegion = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{ Path to an image file. }
  \item{image}{An 'Image' object, if filename is not specified. }
  \item{maxShape}{Maximum size of cell nuclei}
  \item{minShape}{ Minimum size of cell nuclei }
  \item{failureRegion}{ minimum size of failure regions}
}
\details{
  The image is segmented. An image is created, in which every cell is labeled with a numberl. Furthermore, a table including the features of the cells is created.
  In order to create the training set, the table with the cell features has to be opened for instance in a spreadsheet program. Class values for the cells
  have to be inserted in the column 'class'. The corresponding cell in the image can be identified by the column 'index' (numbers in column index correspond to numbers
  in the image). Class values for different classes can be numbers or strings.
  Be careful, this function does not work on MacOsX because of font incompabilities.
}
\value{
  A List containing:
  \item{labeledImage }{Image with labeled cells}
  \item{cellFeatures }{Table of the cell features.}
}
\author{ Henrik Failmezger, failmezger@cip.ifi.lmu.de }


\seealso{ 'createClassifier' }
\examples{
f = system.file("extdata", "exImg.jpg", package="CRImage")
trainingValues=createTrainingSet(filename=f,maxShape=800,minShape=40,failureRegion=2000)
#display(trainingValues[[1]])
#trainingValues[[2]]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
