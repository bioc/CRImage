\name{calculateOtsu}
\alias{calculateOtsu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Does Otsu thresholding
}
\description{
The function applies Otsu thresholding on the image.
}
\usage{
calculateOtsu(allGreyValues)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{allGreyValues}{Vector of grey values.}
}
\details{
The function calculates a value which separates the grey value histogram the best in foreground and background.
}
\value{
the threshold
}
\references{
Nobuyuki Otsu: A threshold selection method from grey level histograms. In: IEEE Transactions on Systems, Man, and Cybernetics. New York 9.1979, S.62-66. ISSN 1083-4419}

\author{
Henrik Failmezger, failmezger@cip.ifi.lmu.de
}


\seealso{
calculateThreshold
localOtsuThreshold

}
\examples{
f1= system.file("extdata", "exImg2.jpg", package="CRImage")
print(f1)
img=readImage(f1)
print(img)
#convert to grayscale
imgG=EBImage::channel(img,'grey')
#threshold value
t=calculateOtsu(as.vector(imgG))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
