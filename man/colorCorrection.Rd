\name{colorCorrection}
\alias{colorCorrection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Color transfer between images.}
\description{
The colors of one image are adapted to the colors of a target image.}
\usage{
colorCorrection(imgO, meanStdTarget,whiteMask = c())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{imgO}{
The image who's colors should be adapted}
  \item{meanStdTarget}{
Array with mean and standard deviation of the target image.}
\item{whiteMask}{Boolean mask of white pixel in the image. These pixels are excluded from color correction.}
}
\details{
Mean and standard deviation of the target image can be calculated using the function calculateMeanStdTarget. }
\value{
The image with adapted colors.
}
\references{
Reinhard, E.; Adhikhmin, M.; Gooch, B.; Shirley, P.; , "Color transfer between images," Computer Graphics and Applications, IEEE , vol.21, no.5, pp.34-41, Sep/Oct 2001
doi: 10.1109/38.946629
URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=946629&isnumber=20481}
\author{
Henrik Failmezger, failmezger@cip.ifi.lmu.de
}



\seealso{
calculateMeanStandardTarget
}
\examples{
#read the target image
	f1= system.file("extdata", "exImg2.jpg", package="CRImage")
	targetImage=readImage(f1)
	#read the image whose color values should be adapted
	f2= system.file("extdata", "exImg3.jpg", package="CRImage")
	imgToConvert=readImage(f2)
	#calculate mean and standard deviation of target color channels
	mst=calculateMeanStdTarget(targetImage)
	# create a white pixel mask
	whitePixelMask=imgToConvert[,,1]>0.85 & imgToConvert[,,2]>0.85 & imgToConvert[,,3]>0.85
	#adapt color channels of image
	imgCorrected=colorCorrection(imgToConvert,mst,whitePixelMask)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
