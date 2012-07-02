\name{labelCells}
\alias{labelCells}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Interactive Session for cell labeling}
\description{
The functions creates an interactive session in order to label cells with their classes.
The labeled cells can be used as training set for the classifier. Note!! This is until now only tested for MacOsX.}
\usage{
labelCells(img, segmentedImage, classes, classColours, nblocks = 3, labeledPoints = NULL, filename = NULL, filenameImage = NULL,transformCoordinates=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{img}{The image.}
  \item{segmentedImage}{The segmented image.}
  \item{classes}{The possible class labels.}
  \item{classColours}{The colors for the class labels.}
  \item{nblocks}{The image can be separated in several blocks, as zooming is not possible.}
  \item{labeledPoints}{Labeled cells from a previous training session.}
  \item{filename}{The table of labeled cells is saved at this location.}
  \item{filenameImage}{The image with the labeled cells is saved at this location.}
\item{transformCoordinates}{deprecated}
}
\details{
Use the keys:
a: In order to add a label to a cell.
d: In order to delete a label from a cell.
c: To switch between classes.
q: To quit the interactive session.
r: To refresh the session (labeled cells will be shown after refreshing)}
\value{
A table with columns:
index: the index of the cell in the segmented image.
x: x-coordinate of the cell
y: y-coordinate of the cell
classCell: Label of the cell
xLocal: Local x coordinate in the subimage(block)
yLocal: Local y coordinate in the subimage(block)
block: Block number in which the cell arises.
}

\author{
Henrik Failmezger}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
\keyword{ misc }% __ONLY ONE__ keyword per line
