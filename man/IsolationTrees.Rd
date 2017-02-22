\name{IsolationTrees}
\alias{IsolationTrees}
\title{Building isolation trees}
\description{
 Building isolation trees
 }
\usage{
       IsolationTrees(x, ntree=10, hlim=as.integer(ceiling(log2(nrow(x)))), rowSamp=F, nRowSamp=nrow(x), nmin=1, rFactor=1, colSamp=F, nColSamp=ncol(x), colWeight=c(rep(1,ncol(x))))
}
\arguments{
  \item{x}{a data frame of training samples}
  \item{ntree}{number of tree to build}
  \item{hlim}{height limit}
  \item{rowSamp}{logical swith to perform random sub-sampling}
  \item{nRowSamp}{sub-sampling size; it must be less than or equal to the training sample size}
  \item{nmin}{ minimum number of sample to form a leaf}
  \item{rFactor}{ randomisation factor, range from 0 to 1, 0 for fully deterministic, 1 for fully random}
  \item{colSamp}{logical switch to perform random attribute-sampling}
  \item{nColSamp}{attribute-sampling size; it must be less than or equal to the number of attributes}
  \item{colWeight}{attribute weight that is being used in random attribute sub-sampling}
}
\value{
a data structure that represent an Isolation Forest model
}
\details{
 Building random binary trees
 }
\references{Fei Tony Liu, Kai Ming Ting, and Zhi-Hua Zhou\cr
 \emph{Isolation Forest}\cr
 IEEE International Conference on Data Mining 2008 (ICDM 08), Pisa, Italy, 2008.
 \url{http://www.gscit.monash.edu.au/gscitweb/loid.php?loid=905282&mimetype=application/pdf}}
\seealso{
  \code{\link{AnomalyScore}}
}
\examples{
library(IsolationForest)
data(stackloss)
# train a model of Isolation Forest
tr<-IsolationTrees(stackloss, rFactor=0)
#evaluate anomaly score
as<-AnomalyScore(stackloss,tr)
# show anomaly score
as$outF
}
\author{Fei Tony Liu}
\keyword{models}


