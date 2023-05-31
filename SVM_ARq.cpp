#include <TMB.hpp>
#include <algorithm>
#include <iostream>

template<class Type>
Type objective_function<Type>::operator() ()
{
  /* Data Section*/
  DATA_VECTOR(Y);
  /* Parameter Section*/
  PARAMETER_VECTOR(X);
  /*PARAMETER_VECTOR(PHI0);*/
  PARAMETER_VECTOR(THETA0);
  PARAMETER(Mu);
  PARAMETER(logSigma);
  PARAMETER(avg);
  
  /* Procedure Section*/  
  int m = X.size();  
  int n = Y.size();
  /*int p = PHI0.size();*/
  int q = THETA0.size();
  int b = m - m;
  Type Sigma = exp(logSigma);
  Type neglogL = 0.0;
  
  /*
  vector<Type> r_phi(p);
  r_phi.fill(0.0);
  
  for(int i = 0; i < p; i++) {
    r_phi(i) = PHI0[i]/(sqrt(1 + pow(PHI0[i],2)));
  }
  
  matrix<Type> y_phi(p,p);
  y_phi.fill(0.0);
  
  for (int k = 0; k<p; k++) {
    for (int i = 0; i < (k); i++) {
      y_phi(k,i) += y_phi(k-1,i) - r_phi(k) * y_phi(k-1,k-i-1);
    }
    y_phi(k,k) += r_phi(k);
  }
  
  vector<Type> PHI = y_phi.row(p-1);
  */
  
  vector<Type> r_theta(q);
  r_theta.fill(0.0);
  
  for(int i = 0; i < q; i++) {
    r_theta(i) = THETA0[i]/(sqrt(1 + pow(THETA0[i],2)));
  }
  matrix<Type> y_theta(q,q);
  y_theta.fill(0.0);
  
  for (int k = 0; k<q; k++) {
    for (int i = 0; i < (k); i++) {
      y_theta(k,i) += y_theta(k-1,i) - r_theta(k) * y_theta(k-1,k-i-1);
    }
    y_theta(k,k) += r_theta(k);
  }
  
  vector<Type> THETA = y_theta.row(q-1);  
  
  
  vector<Type> e(m);
  
  
  for (int i = 0; i < m; i++) {
    Type em = Mu;
    /*
    for (int j = 0; j<p; j++) {
      if (i-1-j < 0) {
        em += 0;
      }
      else {
        em += PHI(j)*(X[i-1-j] - Mu); 
      }
    }
    */
    for (int k = 0; k<q; k++) {
      if (i-1-k < 0) {
        em += 0;
      }
      else {
        em -= THETA(k)*e[i-1-k]; 
      }
      
    }
    neglogL -= dnorm(X[i], em, Sigma, true);
    e(i) = X[i] - em;
  }
  
  
  for (int i = q; i < n; i++) {
    Type e = exp(X[i + b]/2);
    neglogL -= dnorm(Y[i], Y[i-1]+avg, e, true);
  }
  /*REPORT(PHI);*/
  return neglogL;
}
