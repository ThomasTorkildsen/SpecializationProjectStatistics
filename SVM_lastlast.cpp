#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  /* Data Section*/
  DATA_VECTOR(Y);
  /* Parameter Section*/
  PARAMETER_VECTOR(X);
  PARAMETER(Mu);
  PARAMETER(logitPhi);
  PARAMETER(logSigma);
  PARAMETER(Theta)
    /* Procedure Section*/  
    
    int n = Y.size();
  /*vector<Type> X_hat(n);*/
  Type Sigma = exp(logSigma);
  Type Phi = 2 / (1 + exp(-logitPhi)) - 1; 
  Type neglogL = -dnorm(X[0], Mu, Sigma*sqrt(1/(1 - Phi*Phi)),true);
  for (int i = 1; i < n; i++) {
    Type m = Mu + Phi*(X[i-1] - Mu);
    neglogL -= dnorm(X[i], m, Sigma, true);
  }
  for (int i = 1; i < n; i++) {
    Type e = exp(X[i]/2);
    neglogL -= dnorm(Y[i], Y[i-1]+Theta, e, true);
  }
  ADREPORT(Sigma);
  ADREPORT(Phi);
  return neglogL;
}