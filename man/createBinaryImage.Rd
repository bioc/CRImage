\name{createBinaryImage}
\alias{createBinaryImage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Thresholding }
\description{
Creates a binary image from a grayscale image by thresholding.
}
\usage{
createBinaryImage(imgG,img=NULL,method="otsu",threshold=NULL,numWindows=1,whitePixelMask=c())
}
\arguments{
  \item{img}{An Image object or an array.}
\item{imgG}{The grey valued  Image object. }
\item{method}{Either "otsu" or "phansalkar" }
\item{threshold}{Fixed threshold}
\item{numWindows}{Number of windows to use for threshold calculation.}
\item{whitePixelMask}{Boolean mask of white pixels, if they should be excluded from thresholding }
}
\details{
The functions returns the binary image resulting from the thresholding. If threshold is defined, all pixels smaller than this value will be fixed to 1 all other values will be set to 0. 
If threshold is undefined, the thresholding value is calculated automatically using 'otsu' or 'phansalkar' thresholding.

The function 'otsu' does Otsu thresholding on the grey level histograms of the image.
The function 'phansalkar' does thresholding using the mean and standard deviation of a specified window. The thresholding is done on the RGB as well as the LAP color space and the results are ORed.
The window size is dim(img)/numWindows. White pixel can be excluded from thresholding (e.g. if white is background) by defining a whitePixelMask
}
\value{
  The binary image.
}

\author{ Henrik Failmezger, failmezger@lmb.uni-muenchen.de }
\references{
Neerad Phansalkar, Sumit More, Ashish Sabale, Dr. Madhuri Joshi, "Adaptive Local Thresholding for Detection of Nuclei in Diversly Stained Cytology Images," 2011 IEEE International Conference in Communications and Signal Processing (ICCSP), pp. 218, 10 Feb. 2011

Nobuyuki Otsu: A threshold selection method from grey level histograms. In: IEEE Transactions on Systems, Man, and Cybernetics. New York 9.1979, S.62-66. ISSN 1083-4419
}
\examples{
f= system.file("extdata", "exImg.jpg", package="CRImage")
img=readImage(f)
#conversion to grayscale
imgG=EBImage::channel(img,"gray")
imgB=createBinaryImage(imgG,img=img,method="otsu",numWindows=4)
#white pixel mask
whitePixelMask=img[,,1]>0.85 & img[,,2]>0.85 & img[,,3]>0.85
#exclude white pixels from thresholding
imgB=createBinaryImage(imgG,img=img,method="otsu",numWindows=4,whitePixelMask)
#phansalkar threshold
imgB=createBinaryImage(imgG,img=img,method="phansalkar",numWindows=4)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
