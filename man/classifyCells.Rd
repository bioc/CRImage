\name{classifyCells}
\alias{classifyCells}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A function to classify cells }
\description{
  The function classifies cells and paints the different class types in the image.
}
\usage{
classifyCells(classifier,filename="",image=NA,segmentedImage=NA,featuresObjects=NA,paint=TRUE,KS=FALSE,cancerIdentifier=NA, maxShape=NA,minShape=NA,failureRegion=NA,colors=c(),classesToExclude=c(),threshold="otsu",numWindows=2,structures=NA,classifyStructures=FALSE,pixelClassifier=NA,ksToExclude=c())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{classifier}{ A Support Vector Machine created by createClassifier or directly by the package e1071 }
  \item{filename}{ A path to an image file. }
  \item{image}{An 'Image' object or an array.}
  \item{segmentedImage}{An 'Image' object or an array.The corresponding segmented image (created by segmentImage)}
  \item{featuresObjects}{Cell feature file of the segmentedImage (created by segmentImage)}
  \item{paint}{ If true, the classified cells are painted with different colors in the image}
  \item{KS}{ Use Kernel Smoohter in classification? }
  \item{cancerIdentifier}{A string which describes, how the cancer class is named. }
  \item{maxShape}{Maximum size of cell nuclei}
  \item{minShape}{ Minimum size of cell nuclei }
  \item{failureRegion}{ minimum size of failure regions}
  \item{colors}{Colors to paint the classes}
  \item{classesToExclude}{Which class should be excluded?}
  \item{threshold}{Which thresholding method should be used, "otsu" or "phansalkar"}
\item{numWindows}{Number of windows to use for thresholding.}
\item{structures}{If the image is already segmented, structures can be inserted to enable hierarchical classification.}
\item{classifyStructures}{Use hierarchical classification. If yes a pixel classifier has to be defined. }
\item{pixelClassifier}{A SVM to classify pixel based on their color values. Needed if hierarchical classification should be applied.}
\item{ksToExclude}{These classes are excluded from kernel smoothing.}

}
\details{
  The kernels smoother improves the classification for cells which are likely to occur in clusters, like tumour cells. The kernel smoothing method can only be applied for two classes. If there are more classes only the normal svm without kernel smoothing is applied.
Different classes are labeled with different colors in the image.
}
\value{
  A list with
  \item{comp1 }{classes}
  \item{comp2 }{Classes, painted in the image, if paint was true}
}

\author{Henrik Failmezger, failmezger@mpipz.mpg.de}

\examples{
t = system.file("extdata", "trainingData.txt", package="CRImage")
#read training data
trainingData=read.table(t,header=TRUE)
#create classifier
classifier=createClassifier(trainingData)[[1]]
#classify cells
f = system.file("extdata", "exImg.jpg", package="CRImage")
classesValues=classifyCells(classifier,filename=f,KS=TRUE,maxShape=800,minShape=40,failureRegion=2000)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
