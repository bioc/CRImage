\name{convertRGBToHSV}
\alias{convertRGBToHSV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conversion from RGB color space to HSV color space
}
\description{
The RGB Image is converted to an HSV image. }
\usage{
convertRGBToHSV(img)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{img}{The RGB image}
}
\details{
The entries of the array are Hue, Saturation and Value.
}
\value{
The image in HSV color space.
}

\author{
Henrik Failmezger, failmezger@cip.ifi.lmu.de
}


\seealso{
convertHSVToRGB
convertRGBToLAB
convertLABToRGB
}
\examples{
f= system.file("extdata", "exImg.jpg", package="CRImage")
img=readImage(f)
#conversion to HSV color space
imgHSV=convertRGBToHSV(img)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
