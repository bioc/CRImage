\name{correctCopyNumber }
\alias{correctCopyNumber}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Allelic Copy Number correction for cellularity}
\description{
This function segments copy number and corrects log-ratios (LRR) and
beta allele frequencies (BAF) values for cellularity.
}
\usage{
correctCopyNumber(arr="Sample1", chr=NULL, p=NULL, z=NULL, min.value=-5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{arr}{Name of the array.}
  \item{chr}{Chromosome to run. If \code{NULL}, all chromosomes are
    run.}
  \item{p}{Percentage of tumoural cells}.
  \item{z}{Copy Number Data. Must be a dataframe with the following
    columns: Name (id of the probe), Chr (chromosome), Pos (position),
    LRR (log ratios) and BAF (beta allele frequencies).}
  \item{min.value}{Value assigned to the probes that have 0 copies after correction.}
}
\details{
  The data.frame \code{z} must contain only SNP probes, that is probes
  with both LRR and BAF values. It is recommended that all replicated
  probes are merged so the positions are unique.
  This function calls \code{DNAcopy} to segment the LRR and then correct
  the segmented profiles for normal contamination according to the
  method described in the reference below (see for details).
}

\value{
  A list with 2 components:
  \item{y}{a data.frame with as many rows as probes containing the
    following variables: Chrom (chromosome),
    Pos (position), Orig.LRR (LRR before correction) Orig.BAF (BAF
  before correction), Corr.LRR (LRR after cellularity correction) and
  Corr.BAF (BAF after correction)}
  \item{seg}{a data.frame with the segmented data. Contains the
  following columns:    ID (name of the array), chrom (chromosome),
  loc.start (start of the region),  loc.end (end of the region),
  num.mark (number of probes in the region), seg.mean (LRR of the
  region), BAF (BAF of the regions), num.BAF (number of SNP probes in
  the region), Sa (estimated absolute copy number for the first allele),
  Sb (estimated absolute copy number for the first allele), LRR.tum
  (corrected LRR for the region), BAF.tum (corrected BAF for the region).}
}


\note{
Includes an adaptation of \code{aCGH} \code{mergeLevels} function to fix
a problem with \code{ansari.test}.
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


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
