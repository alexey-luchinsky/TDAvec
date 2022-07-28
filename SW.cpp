#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;

// [[Rcpp::export]]
double computeSWdistcpp(NumericMatrix d1, NumericMatrix d2, int homDim, int M=10){
  int n1 = 0; // number of rows with the correct dimension
  for(int i=0;i<d1.nrow();++i){
    if(d1(i,0) == homDim){
      ++n1; 
    }
  }
  
  int n2 = 0; // number of rows with the correct dimension
  for(int i=0;i<d2.nrow();++i){
    if (d2(i,0) == homDim){
      ++n2; 
    }
  }
  
  NumericVector x1(n1),y1(n1);
  int n=0;
  for(int i=0;i<d1.nrow();++i){
    if (d1(i,0) == homDim){
      x1[n] = d1(i,1);
      y1[n] = d1(i,2);
      ++n;
    }
  }
  NumericVector p1=(x1+y1)/2;
  
  NumericVector x2(n2),y2(n1);
  n=0;
  for(int i=0;i<d2.nrow();++i){
    if( d2(i,0) == homDim){
      x2[n] = d2(i,1);
      y2[n] = d2(i,2);
      ++n;
    }
  }
  NumericVector p2=(x2+y2)/2;
  
  double theta = -M_PI_2;
  double s = M_PI/M;
  
  NumericMatrix V1(n1+n2,M);
  NumericMatrix V2(n1+n2,M);
  
  for(int m=0;m<M;++m){
    for(int i=0;i<(n1+n2);++i){
      if (i<n1){
        V1(i,m) = x1[i]*cos(theta) + y1[i]*sin(theta);
      } else{
        V1(i,m) = p2[i-n1]*(cos(theta) + sin(theta));
      }
      if (i<n2){
        V2(i,m) = x2[i]*cos(theta) + y2[i]*sin(theta);
      } else{
        V2(i,m) = p1[i-n2]*(cos(theta) + sin(theta));
      } 
    }
    theta = theta + s;
  }
  
  NumericVector v1(n1+n2);
  NumericVector v2(n1+n2);
  
  double SW=0;
  for(int m=0;m<M;++m){ 
    v1 = V1(_,m);
    v2 = V2(_,m);
    SW = SW + sum(abs(v1.sort()-v2.sort()))/M;
  }
return SW;
}

