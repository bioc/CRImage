\name{plotCorrectedCN}
\alias{plotCorrectedCN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot CN profiles corrected for cellularity}
\description{
This function takes the result of a call to \code{correctCopyNumber} and
plots the results.
}
\usage{
plotCorrectedCN(CN, chr=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CN}{object result of a call to \code{correctCopyNumber}.}
  \item{chr}{chromosome to plot.}
}
\details{
  A panel with four plots is created. The top panel shows LRR (with
  DNAcopy segmentation overlayed) and BAF before correction and the
  bottom panel shows the plots after correction.
}

\value{
No value is returned.
}


\references{
Yuan, Y et al. Quantitative image analysis of cellular heterogeneity in
primary breast tumors enriches genomic assays. In prep.
}

\author{Oscar M. Rueda, rueda.om@gmail.com}

\examples{
LRR <- c(rnorm(100, 0, 1), rnorm(10, -2, 1), rnorm(20, 3, 1),
       rnorm(100,0, 1)) 
BAF <- c(rnorm(100, 0.5, 0.1), rnorm(5, 0.2, 0.01), rnorm(5, 0.8, 0.01), rnorm(10, 0.25, 0.1), rnorm(10, 0.75, 0.1),
       rnorm(100,0.5, 0.1)) 

Pos <- sample(x=1:500, size=230, replace=TRUE)
Pos <- cumsum(Pos)
Chrom <- rep(1, length(LRR))
z <- data.frame(Name=1:length(LRR), Chrom=Chrom, Pos=Pos, LRR=LRR, BAF=BAF)
res <- correctCopyNumber(arr="Sample1", chr=1, p=0.75, z=z)
plotCorrectedCN(res, chr=1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
