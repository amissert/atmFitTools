#ifndef T2KSAMPLE_H 
#define T2KSAMPLE_H 

#include "TH2FV.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TString.h"
#include "TFile.h"
#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TText.h"
#include "TLine.h"
#include <iostream>

#define NMCMCMAX 200
#define TOTALBINSMAX 50

const double EnuBinningMuons[] = {0.,500,900,1500,30000};
const int NEnuBinsMuons = 4;
const double EnuBinningElectrons[] = {0.,400,800,1200.};
const int NEnuBinsElectrons = 3;

class t2kSample{
  public:

  t2kSample(const char* nametag, int flgelectron, int flgmuon);
  TString nameTag;

  TH1D* hEnuElectron;
  TH1D* hEnuMuon;
  float enuBinValues[NMCMCMAX][TOTALBINSMAX];
  float enuBinMean[TOTALBINSMAX];
  float enuBinCov[TOTALBINSMAX];
  float enuBinCor[TOTALBINSMAX];
  TH2D* corMatrix;
  TH2D* covMatrix;
  TH1D* binMean;
  int toyCounter;
  int nTotBins;
  int flgE;
  int flgMu;

  void resetAll();
  void initHistos();
  void finishToyRun(); 
  void saveHistos(const char* directory);
  void readHistos();
  void calcUncertainties();

};

#ifdef CINTMODE
#include "t2kSample.cxx"
#endif

#endif
