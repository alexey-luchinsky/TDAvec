#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;

//' Calculates the Persistence Landscape
//' 
//' @param D N by 3 matrix (columns contain dimension, birth and persistence values respectively)
//' @param homDim homological dimension (0 for H0, 1 for H1, etc.)
//' @param k order of landscape function
//' @param scaleSeq sequence of scale values for vectorization
//' @examples
//' N <- 100
//' set.seed(123)
//' X <- TDA::circleUnif(N) + rnorm(2*N,mean = 0,sd = 0.2)
//' D <- TDA::ripsDiag(X,maxdimension = 1,maxscale = 2)$diagram 
//' scaleSeq = seq(0,2,length.out=11) # sequence of scale values
//' computePL(D,homDim=0,k=1,scaleSeq)
// [[Rcpp::export]]
NumericVector computePL(NumericMatrix D, int homDim, int k, NumericVector scaleSeq){
  int n_rows = 0; // number of rows with the correct dimension
  for(int i=0;i<D.nrow();++i){
    if(D(i,0) == homDim){
      ++n_rows; 
    }
  }
  
  int L = scaleSeq.size();
  if (n_rows == 0) return NumericVector(L);
  
  NumericVector x(n_rows),y(n_rows);
  int n=0;
  for(int i=0;i<D.nrow();++i){
    if( D(i,0) == homDim){
      x[n] = D(i,1);
      y[n] = D(i,2);
      ++n;
    }
  }
  
  NumericVector lambda(L);
  for (int i=0;i<L;++i){
    NumericVector Lambda = pmax(pmin(scaleSeq[i] - x, y - scaleSeq[i]),0);
    lambda[i] = Lambda.sort(true)[k-1];
  }
  return lambda; 
}

