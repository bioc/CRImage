\name{convertRGBToLAB}
\alias{convertRGBToLAB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Converts RGB to LAB color space.
}
\description{
Conversion of Color spaces.
}
\usage{
convertRGBToLAB(imgT)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{imgT}{The RGB image.}
}
\details{
Color space conversion}
\value{
The image in LAB color space.}

\author{
Henrik Failmezger, failmezger@cip.ifi.lmu.de }



\examples{
f= system.file("extdata", "exImg.jpg", package="CRImage")
img=readImage(f)
#conversion to LAB color space
imgLAB=convertRGBToLAB(img)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
\keyword{ misc}% __ONLY ONE__ keyword per line
