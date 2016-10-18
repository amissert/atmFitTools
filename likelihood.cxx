#ifndef LIKELIHOOD_CXX
#define LIKELIHOOD_CXX

#include "likelihood.h"

using namespace std;

/////////////////////////////////////////
//likelihood from tn186
double evalLnL(double ndata, double nmc, double norm){
  if (nmc==0) return 0.;
  if (ndata==0) return nmc*norm;
  double Nmc = nmc*norm;
  double ngLnL = (Nmc-ndata) + ndata*TMath::Log(ndata/Nmc);
  return ngLnL;
}

/////////////////////////////////////////////////////////////////////////
// Gaussian with MC errors
double evalGausChi2WithError(double ndata, double mcmean, double mcsig){
  double sigmasq = ndata + (mcsig*mcsig);
  double diff = (ndata-mcmean);
  return (diff*diff)/(2.*sigmasq);

}


////////////////////////////////////////////////////////////////
// standard chi2
double evalChi2(double ndata, double mcexpect){
  double diff = (ndata-mcexpect);
  return (diff*diff)/(2.*ndata);
}

#endif
