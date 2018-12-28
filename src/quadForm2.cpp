#include <Rcpp.h>
#include <Rmath.h>
#include <R.h>

using namespace Rcpp;


// [[Rcpp::export]]
double xQx1(NumericMatrix Q, NumericVector x){
int n = Q.nrow();  
double phiQphi = 0;
double Qphi = 0;
for(int i = 0; i < n; i++) {
  Qphi = 0;
 for(int j = 0; j < n; j++) {
      Qphi +=  Q(i, j) * x[j];
    }
    phiQphi += Qphi*x[i];
}
return(phiQphi);
}