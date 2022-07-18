#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;
using namespace std;

NumericVector seq_C(double a, double b, double by){
  int n = floor((b-a)/by)+1;
  NumericVector out(n);
  for(int i=0; i<n; ++i){
    out[i] = a + i*by;
  }
  return out;
}

NumericVector rowSums_C(NumericMatrix D) {
  NumericVector out(D.nrow());
  for ( int i = 0; i < D.nrow(); ++i ) {
    out[i] = sum(D(i,_));
  }
  return out;
}

NumericVector outer_C(NumericVector x, NumericVector y) {
  int n = x.size();
  int m = y.size();
  NumericVector result(n*m);
  for ( int i = 0; i < m; ++i ) {
    for ( int j = 0; j < n; ++j ) {
      result[j+i*n] = x[j] * y[i];
    }
  }
  return result;
}

// PS for H0
NumericVector PSurfaceH0(NumericVector point,
                         NumericVector y_lower,
                         NumericVector y_upper,
                         double sigma,
                         double maxP) {
  double y = point[1];
  NumericVector out2 = pnorm(y_upper,y,sigma,true,false) - pnorm(y_lower,y,sigma,true,false);
  double wgt = y/maxP * (y < maxP) + 1 * (y >= maxP);
  return wgt*out2;
}

// PS for Hk, k>0
NumericVector PSurfaceHk(NumericVector point,
                         NumericVector y_lower,
                         NumericVector y_upper,
                         NumericVector x_lower,
                         NumericVector x_upper,
                         double sigma,
                         double maxP) {
  double x = point[0];
  double y = point[1];
  NumericVector out1 = pnorm(x_upper,x,sigma,true,false) - pnorm(x_lower,x,sigma,true,false);
  NumericVector out2 = pnorm(y_upper,y,sigma,true,false) - pnorm(y_lower,y,sigma,true,false);
  double wgt = y/maxP * (y < maxP) + 1 * (y >= maxP);
  return wgt*outer_C(out1,out2); 
}

// [[Rcpp::export]]
NumericVector computePIcpp(NumericMatrix D,int homDim,
                           int res, double sigma,
                           double minB, double maxB,
                           double minP, double maxP){
// D - N by 3 matrix (columns contain dimension, birth and persistence values respectively)
int n_rows = 0; // number of rows with the correct dimension
  for(int i=0; i<D.nrow(); ++i) {
    if(D(i,0) == homDim) {
      n_rows = n_rows+1; 
    }
  }
  
if (n_rows == 0) return NumericVector(res*res);

NumericMatrix D_(n_rows,2);
int j=0;
for(int i=0; i<D.nrow(); ++i) {
  if(D(i,0) == homDim) {
    D_(j,0) = D(i,1);
    D_(j,1) = D(i,2);
    ++j;
  }
}

double dy = (maxP-minP)/res;
NumericVector y_lower = seq_C(minP,maxP-dy,dy);
NumericVector y_upper = y_lower + dy;
int size;
if ((sum(abs(D_(_,0)))==0)&(homDim==0)) {
  size = res;
  }else{
  size = res*res;
}
NumericMatrix Psurf_mat(size,n_rows);
  
if (size==res) {
  for(int i=0; i<n_rows; ++i){
    Psurf_mat(_,i) = PSurfaceH0(D_(i,_),y_lower,y_upper,sigma,maxP);
  }

} else{
  double dx = (maxB-minB)/res;
  NumericVector x_lower = seq_C(minB,maxB-dx,dx);
  NumericVector x_upper = x_lower + dx;
  for(int i=0; i<n_rows; ++i){
    Psurf_mat(_,i) = PSurfaceHk(D_(i,_),y_lower,y_upper,x_lower,x_upper,sigma,maxP);
  }
}
NumericVector out = rowSums_C(Psurf_mat);
return out;
}