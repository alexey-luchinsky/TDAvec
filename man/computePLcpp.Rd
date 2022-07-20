\name{computePLcpp}
\alias{computePLcpp}
\title{Persistense Landscape}
\usage{
  computePLcpp(D, homDim, k, scaleSeq)
}
\arguments{
  \item{D}{N by 3 matrix (columns contain dimension, birth and persistence values respectively)}
  \item{homDim}{homological dimension (0 for H0, 1 for H1, etc.)}
  \item{k}{order of landscape function}
  \item{scaleSeq}{sequence of scale values for vectorization}
}
\description{
  Calculates the Persistence Landscape, C++ version
}
\examples{
N <- 100
set.seed(123)
X <- TDA::circleUnif(N) + rnorm(2*N,mean = 0,sd = 0.2)
D <- TDA::ripsDiag(X,maxdimension = 1,maxscale = 2)$diagram 
scaleSeq = seq(0,2,length.out=11) # sequence of scale values
computePLcpp(D,homDim=0,k=1,scaleSeq)
}
