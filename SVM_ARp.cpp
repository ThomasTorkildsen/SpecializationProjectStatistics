#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  /* Data Section*/
  DATA_VECTOR(Y);
  /* Parameter Section*/
  PARAMETER_VECTOR(X);
  PARAMETER_VECTOR(PHI0);
  /*PARAMETER_VECTOR(THETA);*/
  PARAMETER(Mu);
  PARAMETER(logSigma);
  PARAMETER(avg);
  
  /* Procedure Section*/  
  int m = X.size();  
  int n = Y.size();
  int p = PHI0.size();
  int b = m - n;
  /*int q = THETA.size();*/
  Type Sigma = exp(logSigma);
  Type neglogL = 0.0;
  
  vector<Type> r(p);
  r.fill(0.0);
  /*for(int i = 0; i < p; i++) {
   r[i] = 2/(1 + exp(-PHI0[i])) - 1;
  }*/
  
  for(int i = 0; i < p; i++) {
    r(i) = PHI0[i]/(sqrt(1 + pow(PHI0[i],2)));
  }
  matrix<Type> y(p,p);
  y.fill(0.0);
  
  for (int k = 0; k<p; k++) {
    for (int i = 0; i < (k); i++) {
      y(k,i) += y(k-1,i) - r(k) * y(k-1,k-i-1);
    }
    y(k,k) += r(k);
  }
  
  vector<Type> PHI = y.row(p-1);
  
  /*for(int j = 0; j < p; j++) {
   X[j] = Mu;
  }*/
  
  
  for (int i = 0; i < m; i++) {
    Type em = Mu;
    for (int j = 0; j<p; j++) {
      if (i-1-j < 0) {
        em += 0;
      }
      else {
        em += PHI(j)*(X[i-1-j] - Mu); 
      }
    }
    /*    for (int k = 0; k<q; k++) {
     m += THETA[k]*e[i-1-k];
    }*/
    neglogL -= dnorm(X[i], em, Sigma, true);
  }
  
  
  for (int i = p; i < n; i++) {
    Type e = exp(X[i + b]/2);
    neglogL -= dnorm(Y[i], Y[i-1]+avg, e, true);
  }
  /*REPORT(PHI);*/
  return neglogL;
}
