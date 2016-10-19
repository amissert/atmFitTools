#ifndef CALCRES_H
#define CALCRES_H

#include "TH1D.h"

#include <iostream>
#include "TMath.h"

using namespace std;


void cumulatebackward(TH1D* hh);
double integrateError(TH1D* hh, int binlow, int binhigh);
void cumulateforward(TH1D* hh);
int getMostSigBin(TH1D* hh, double baseline=0.0);
void calcCumDiff(TH1D* hmc, TH1D* hdata);
void calcCumFracDiff(TH1D* hmc, TH1D* hdata);
double calcFracDiff(TH1D* hmc, TH1D* hdata);
double calcCumFracDiffI(TH1D* hmc, TH1D* hdata);
double calcResErrorI(TH1D* hmc, TH1D* hdata);
double calcResErrorFrac(TH1D* hmc, TH1D* hdata);
void absVal(TH1D* hh);
double calcResErrorFrac2(TH1D* hmc, TH1D* hdata);
double calcResErrorFrac3(TH1D* hmc, TH1D* hdata);
void getSignif(TH1D* hmc, TH1D* hdata);
double calcResError(TH1D* hmc, TH1D* hdata);
void applyResError(TH1D* hmc, TH1D* hdata);

#ifdef CINTMODE
#include "calcResErr.cxx"
#endif


#endif




