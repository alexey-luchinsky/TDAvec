\name{computePIcpp}
\alias{computePIcpp}
\title{Persistense Image, C version}
\usage{
  computePIcpp(D,homDim,res,sigma,minB,maxB,minP,maxP)
}
\arguments{
  \item{D}{this is a N by 3 matrix (columns contain dimension, birth and persistence values respectively)}
  \item{homDim}{homological dimension (0 for H0, 1 for H1, etc.)}
  \item{res}{resolution parameter}
  \item{sigma}{sigma parameter}
  \item{minB}{minimal birth value}
  \item{maxB}{maximal birth value}
  \item{minP}{minimal persistance value}
  \item{maxP}{maximal persistance value}
}
\description{
  Calculates the Persistence Image, C version
}
\examples{
  N <- 100
  set.seed(123)
  X <- TDA::circleUnif(N) + rnorm(2*N,mean = 0,sd = 0.2)
  # compute PD using Rips filtration
  D <- TDA::ripsDiag(X,maxdimension = 1,maxscale = 2)$diagram
  D[,3] <- D[,3] - D[,2]
  colnames(D)[3] <- "Persistence"
  res <- 5 # resolution or grid size
  minPH0 <- min(D[D[,1]==0,3]); maxPH0 <- max(D[D[,1]==0,3])
  sigma <- 0.5*(maxPH0-minPH0)/res # default way of selecting sigma - can be overridden
  computePIcpp(D,homDim=0,res,sigma,minB=NA,maxB=NA,minPH0,maxPH0)
}
