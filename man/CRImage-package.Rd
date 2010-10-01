\name{CRImage-package}
\alias{CRImage-package}
\alias{CRImage}
\docType{package}
\title{
CRImage is a package to analyze images and classify cells.
}
\description{
CRImage allows classification of cells in biological images. It offers methods to segment cells or cell nuclei
in biological images for example HE stained images. It offers methods to create a classifier and to classify cells in these images.
Furthermore it allows the calculation of tumour cellularity for large microscope images.

CRImage makes use of the image processing package EBImage, which uses the 'ImageMagick' library for image I/O operations and the 'GTK' library to display images.
}
\details{
\tabular{ll}{
Package: \tab CRImage\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-04-27\cr
License: \tab  LGPL Version 2 or later\cr
LazyLoad: \tab yes\cr
}
}
\section{Package content}{

Image processing methods:
\itemize{
\item{calculateThreshold}
\item{segmentImage}
}
Classification:
\itemize{
\item{createTrainingSet}
\item{createClassifier}
\item{classifyCells}
}
Tumour cellularity
\itemize{
\item{calculateCellularity}
\item{processAperio}
}
}
\author{
Henrik Failmezger, <failmezger@cip.ifi.lmu.de>

Yinyin Yuan, <Yinyin.Yuan@cancer.org.uk>

Oscar Rueda, <oscar.rueda@cancer.org.uk>

Florian Markowetz, <florian.markowetz@cancer.org.uk>

CRI Cambridge

Li Ka Shing Centre

Robinson Way

Cambridge, CB2 0RE, UK

Ludwigs-Maximilians University of Munich


}


\keyword{misc}

\examples{
example(segmentImage)
example(createClassifier)
example(classifyImage)
}
