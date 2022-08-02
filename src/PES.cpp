#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;

//' computePES
//'
//' @examples
//' N <- 100
//' set.seed(123)
//' X <- TDA::circleUnif(N) + rnorm(2*N,mean = 0,sd = 0.2)  
//' D <- TDA::ripsDiag(X,maxdimension = 1,maxscale = 2)$diagram  # compute PD using Rips filtration
//' scaleSeq = seq(0,2,length.out=11) # sequence of scale values
//' computePES(D,homDim = 0,scaleSeq) # compute PES for homological dimension H0
//' computePES(D,homDim = 1,scaleSeq) # compute PES for homological dimension H1
// [[Rcpp::export]]
NumericVector computePES(NumericMatrix D, int homDim, NumericVector scaleSeq){
  int n_rows = 0; // number of rows with the correct dimension
  int scaleLen = scaleSeq.size()-1;
  for(int i=0;i<D.nrow();++i){
    if((D(i,0) == homDim)&(Rcpp::traits::is_finite<REALSXP>(D(i,2)))){
      ++n_rows; 
    }
  }
  
  if (n_rows == 0) return NumericVector(scaleLen);
  
  NumericVector x(n_rows),y(n_rows);
  int n=0;
  for(int i=0;i<D.nrow();++i){
    if((D(i,0) == homDim)&(Rcpp::traits::is_finite<REALSXP>(D(i,2)))){
      x[n] = D(i,1);
      y[n] = D(i,2);
      ++n;
    }
  }
 
  NumericVector lL = (y - x)/sum(y-x);
  NumericVector entr = -lL*log10(lL)/log10(2);
  
  NumericVector pes(scaleLen);
  NumericVector b(n);
  for (int k=0;k<scaleLen;++k){
    b = pmin(scaleSeq[k+1],y)-pmax(scaleSeq[k],x);
    pes[k] = sum(entr*pmax(0,b))/(scaleSeq[k+1]-scaleSeq[k]);
  }
  return pes; 
}

