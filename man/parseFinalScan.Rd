\name{parseFinalScan}
\alias{parseFinalScan}
\title{ Internal function to parse CWS index files.
}
\description{
The function parses the FinalScan.ini file, to a R readable format.}
\usage{
parseFinalScan(file)
}
\arguments{
  \item{file}{
	CWS file
}
}
\details{
The FinalScan.ini files keep information about the information of subimages in the whole image.}
\value{
 List with subimage names and corresponding slides.	
}

\author{
Henrik Failmezger, <email: failmezger@cip.ifi.lmu.de>
}
\examples{
t = system.file("extdata", "trainingData.txt", package="CRImage")
#read training data
trainingData=read.table(t,header=TRUE)
#create classifier
classifier=createClassifier(trainingData,topo=FALSE)[[1]]
#classify aperio
f = system.file("extdata",  package="CRImage")
f=file.path(f,"8905")
dir.create("AperiOutput")
#takes long time!
#processAperio(classifier=classifier,inputFolder=f,outputFolder="AperiOutput",identifier="Da",numSlides=2,cancerIdentifier="1",maxShape=800,minShape=40,failureRegion=2000)

}
\keyword{ misc }
