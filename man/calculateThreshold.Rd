\name{calculateThreshold}
\alias{calculateThreshold}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Thresholding }
\description{
  Calculates the grey value which separates the grey-level image best in foreground and background.
  Tthe Otsu Method is used for calculating the threshold.
}
\usage{
calculateThreshold(allGreyValues)
}
\arguments{
  \item{allGreyValues}{a vector of grey values. }
}
\details{
The optimal threshold is searched by a histogram separation method.
}
\value{
  The calculated threshold.
}
\references{ Otsu, N. A threshold selection method from gray level histograms IEEE Trans. Systems, Man and Cybernetics, 1979, 9, 62-66}
\author{ Henrik Failmezger, failmezger@cip.ifi.lmu.de }

\examples{
f= system.file("extdata", "exImg.jpg", package="CRImage")
img=readImage(f)
#convert to grayscale
imgB=channel(img,"gray")
#find white pixels and exclude them from thresholding(if white is background)
indexWhitePixel=which(img[,,1]>0.85 &img[,,2]>0.85 & img[,,3]>0.85)
#calculate threshold
t=calculateThreshold(as.vector(imgB[-indexWhitePixel]))
#create binary image
imgB[imgB>t]=-1
imgB[imgB != -1]=0
imgB[imgB == -1]=1
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
