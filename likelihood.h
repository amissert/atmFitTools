#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H
#include "TMath.h"

#include "TH1D.h"
#include "TGraph.h"

#include <iostream>
#include <math.h>

double evalLnL(double ndata, double nmc, double norm = 1.);
double evalChi2(double ndata, double mcexpect);
double evalGausChi2WithError(double ndata, double mcmean, double mcsig);

#ifdef CINTMODE
#include "likelihood.cxx"
#endif

#endif
