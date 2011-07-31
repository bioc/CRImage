\name{convertLABToRGB}
\alias{convertLABToRGB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conversion of LAB colour space to RGB colour space
}
\description{
Color space conversion.}
\usage{
convertLABToRGB(imgLAB)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{imgLAB}{
LAB channel vectors.
}
}
\details{
Color space conversion}
\value{
RGB channel vectors.
}

\author{
Henrik Failmezger, failmezger@cip.ifi.lmu.de
}




\examples{
f= system.file("extdata", "exImg.jpg", package="CRImage")
img=readImage(f)
#conversion to HSV color space
imgRGB=convertLABToRGB(img)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc}
\keyword{ misc }% __ONLY ONE__ keyword per line
