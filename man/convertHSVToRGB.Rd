\name{convertHSVToRGB}
\alias{convertHSVToRGB}
\title{
Conversion from HSV color space to RGB color space
}
\description{
The function converts images in the HSV colour space to the RGB colour space.
}
\usage{
convertHSVToRGB(imgHSV)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{imgHSV}{An 'Image' object or an array in the HSV colour space.}
}
\details{
Standard colour space conversion.
}
\value{
An array in the RGB colour space.
}

\author{
Henrik Failmezger, failmezger@cip.ifi.lmu.de
}



\seealso{
convertRGBToHSV
convertRGBToLAB
convertLABToRGB
}
\examples{
f= system.file("extdata", "exImg.jpg", package="CRImage")
img=readImage(f)
#conversion to RGB color space
imgRGB=convertHSVToRGB(img)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
