\name{computeVPBcpp}
\alias{computeVPBcpp}
\title{Vector Persistense Blocks, C version}
\usage{
computeVPBcpp(D, homDim, xSeq, ySeq, tau)
}
\arguments{
\item{D}{N by 3 matrix (columns contain dimension, birth and persistence values respectively)}
\item{homDim}{homological dimension (0 for H0, 1 for H1, etc.)}
\item{xSeq}{sequence of x (birth) values of the grid vertices}
\item{ySeq}{sequence of y (persistence) values of the grid vertices}
\item{tau}{parameter (between 0 and 1) controlling block width. weight function w(x,y) = x+y}
}
\description{
Calculates the Vector Persistence Block
}
\examples{
  X <- TDAstats::circle2d + rnorm(200,mean = 0,sd = 0.1)
  D <- TDAstats::calculate_homology(X,dim = 1,threshold = 2)
  D[,3] <- D[,3] - D[,2]
  xSeq <- unique(quantile(D[D[,1]==1,2],probs = seq(0,1,by=0.1)))
  ySeq <- unique(quantile(D[D[,1]==1,3],probs = seq(0,1,by=0.1)))
  computeVPBcpp(D,homDim = 0, xSeq=xSeq, ySeq,tau = 0.5)
}
