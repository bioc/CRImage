\name{createClassifier}
\alias{createClassifier}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Construction of a classifier }
\description{
  Creates a classifier for a training set.
}
\usage{
createClassifier(trainingData, cross = FALSE, topo = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trainingData}{ A table, created by segmentImage with manually added classes. }
  \item{cross}{ Does 10-fold cross validation to test the classifiers performance.}
  \item{topo}{ Use topological features.}
}
\details{
  Topological features include the density of cells and the size of the surrounding cytoplasma of a cell. These features depend on the size of the image. If training image and the image to classify have different size, these features can fool the classification and should not be enabled.
}
\value{
  A List containing:
  \item{classifier }{The classifier}
  \item{performance }{cross validation performance}
}
\author{ Henrik Failmezger, failmezger@cip.ifi.lmu.de }

\seealso{ 'createTrainingSet','classifyCells' }
\examples{
f = system.file("extdata", "trainingData.txt", package="CRImage")
#read training data
trainingData=read.table(f,header=TRUE)
#create classifier
classifier=createClassifier(trainingData,topo=FALSE)[[1]]

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
