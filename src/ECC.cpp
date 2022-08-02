#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;

NumericVector rowSums_C(NumericMatrix D);
NumericVector computeVAB(NumericMatrix D, int homDim, NumericVector scaleSeq);

//' computeECC
//' 
//' @examples
//' N <- 100
//' set.seed(123)
//' X <- TDA::circleUnif(N) + rnorm(2*N,mean = 0,sd = 0.2)
//' D <- TDA::ripsDiag(X,maxdimension = 1,maxscale = 2)$diagram 
//' scaleSeq = seq(0,2,length.out=11)
//' computeECC(D,maxhomDim = 1,scaleSeq)
// [[Rcpp::export]]
NumericVector computeECC(NumericMatrix D, int maxhomDim, NumericVector scaleSeq){
  NumericMatrix ecc(scaleSeq.size()-1,maxhomDim+1);
  for (int d=0;d<=maxhomDim;++d){
    ecc(_,d) = pow(-1,d)*computeVAB(D,d,scaleSeq);
  }
  return rowSums(ecc);
}

