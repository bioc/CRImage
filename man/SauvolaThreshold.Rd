\name{SauvolaThreshold}
\alias{SauvolaThreshold}
\title{
Do Sauvola thresholding
}
\description{
Thresholding method using mean and standard deviation.
}
\usage{
SauvolaThreshold(allGreyValues)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{allGreyValues}{Vector of gray values.}
}
\details{
A threshold for the gray values is returned
}
\value{
The threshold.
}
\references{
J. Sauvola, M. Pietikainen, "Adaptive Document Image Binarization," Pattern Recognition, vol. 33, 225-236, 2000
}
\author{
Henrik Failmezger, failmezger@cip.ifi.lmu.de
}



\seealso{
createBinaryImage
}
\examples{
f1= system.file("extdata", "exImg2.jpg", package="CRImage")
print(f1)
img=readImage(f1)
print(img)
#convert to grayscale
imgG=EBImage::channel(img,'grey')
#threshold value
t=SauvolaThreshold(as.vector(imgG))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
