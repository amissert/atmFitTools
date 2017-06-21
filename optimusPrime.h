#ifndef OPTIMUS_H
#define OPTIMUS_H

#define CINTMODE 

#include "fqProcessedEvent.h"
#include "mcLargeArray.h"
#include "randomList.h"
#include "TH2FV.h"
#include "TMath.h"
#include "TFile.h"
#include "TChain.h"
#include "TGraph2D.h"
#include "TCanvas.h"
#include "eventSelectors.h"
#include "moreUncertainties.h"
#include "mcmcApply.h"
#include <exception>
#include "TTree.h"
#include "THStack.h"
#include "TRandom2.h"
#include "TPaveText.h"
#include "defineSmall.C"
#include "atmFitPars.h"
#include "stats.C"
#include "smartPlot.C"
#include "summaryPlots.h"
#include <iostream>
#include <algorithm>

#define PRINTSUMMARY
#define NFOMBINSMAX 100


double EnuBinningMuonTh23[] = {0.,2000.};
int    EnuNBinsMuonTh23 = 1;

double EnuBinningMuondM2[] = {0.,600.,2000.};
int    EnuNBinsMuondM2 = 2;

double EnuBinningElectron[] = {0.,1250.};
int    EnuNBinsElectron = 1;


////////////////////////////////////////
// Class to choose FV optimal point
class optimusPrime{

  public:

  // constructor 
  optimusPrime(TChain* t2kmc,
               int nevents, 
               const char* datadir,
               const char* nuefile,
               const char* numufil,
               const char* nue1rpifile);
 
  // vars
  TChain* chmc; //< T2K MC events
  fqProcessedEvent* mcevent; //< individual event values
  int nevents;  //< total number of events to use
  randomList* eventlist; //< random sampling of events if we don't want to use all
  mcLargeArray* fastevents; //< very large array for all T2K FC events
  moreUncertainties* uncertaintyCalculator; //< can calculate addtional uncertainties from entering bg, ect.

  // separate for each:
  TString UncMapFileNameNuE;
  TString UncMapFileNameNuMu;
  TString UncMapFileNameNuE1RPi;
  moreUncertainties* uncNuE;
  moreUncertainties* uncNuE1RPi;
  moreUncertainties* uncNuMu;

  mcmcApply* modifier; //< applies atm fit pars to a given event
  fqcutparams cutPars; //< cut parameters structure defined in eventSelectors.h
//  TString mcmcParFileName;
  TString cardFileName;
  TString pltTag;
  TString outDir;
  summaryPlots* plots;
  summaryPlots* plots1;
  summaryPlots* plots2;
  TPaveText* txtSummary;

//  TRandom2* randy();

  TH2FV* hFV;
  TH2FV* hFVMask;
  TH2FV* hFVMaskSg;
  TH2FV* hFVMaskBg;
  TH2FV* hDeltaSg;
  TH2FV* hDeltaBg;
  TH2FV* hDelta;
  TGraph2D* grDelta;
  TH2FV* hFVAvg;
  TH2FV* hFVAll;
  TH1D*  hErecSeedNuMudM2;
  TH1D*  hErecSeedNuMuTh23;
  TH1D*  hErecSeedNuE;
  TH1D*  hErec[10];
  TH2FV* hFVSummary[8];
  TH1D*  hSummary[20];
  TH1D*  hCurve;
  TCanvas* multiPad;
  TCanvas* canPad;
  THStack* hs;

  // set these to the appropriate attribute[] index
  int indexPIDPar;
  int indexPi0Par;
  int indexPiPPar;
  int indexRCPar;

  int flgPrintSummary;

  float DeltaSg;
  float DeltaBg;
  float Delta;

  float Nevents;
  float Power;
  float Syst;
  float NS;
  float NB;
  float Scale;
  float SysScale;
  int FOMType;
  int AvgType;
  int flgUseSpectrum;
  int bestFOMbin;
  int flgUseFloor;
  int flgUseCCnQE;
  float bestFOMvalue;
  // methods
 
  // masks
  int flgUseMask;
  float maskThresh;
  TH2FV* hMask;
  TH2FV* hMaskOne;
  float smallVariation; 


  void fillFVHistoFast(); 

  float calcFOMBinned(int nselection, float twcut, float wcut, int oscpar, int iplt = 1);

  float calcFOM(float* pow, float* nev, float* sys, int nbin);
  void calcFOMMap(float towallmax, float wallmax,int oscpar, int npts, int nselection);
  float getOscPower(int nutype, int oscpar);
  float getEventWeight(int iev);
  float getOscPowerFast(int nutype, int ievent, int oscpar);
  float getSystUncertainty(int iev,int nutype=14);
  void fillArray();
  void calcFVSummary(int oscpar, int nutype=14);
  int passNuMuCuts(int iev);
  int passNuECuts(int iev);
  double calcDeltaNuMu(float tw1, float w1, float tw2, float w2);
  void calcDeltaMapNuMu(float twbest, float wbest, float twmax, float wmax, int npts=15);
  double calcDeltaNuE(float tw1, float w1, float tw2, float w2);
  void calcDeltaMapNuE(float twbest, float wbest, float twmax, float wmax, int npts=15);

  int isSmallDifference(float tw1, float w1, float tw2, float w2, int oscpar, int nselection);

  // superseeds previous methods for applying numu or nue cuts to a modified event
  int applyCutsToModifiedEvent(int iev);
  /////////////////////////////////////////////////////////////////////////////////
 
  void makeAllPlots(float twmax, float wmax, int oscpar, int npts,int nselection);

  void compareCuts(float tw1, float w1, float tw2, float w2, int oscpar, int nselection);
  void compareFOM(float tw1, float w1, float tw2, float w2, int oscpar, int nselection);
  void showBreakdown();
  void printCutDiff(int nselection);
  void printCompare(const char* dir,float tw1, float w1, float tw2, float w2, int oscpar, int nselection);

  void deleteHistos();
  void initHistos(int nselection, int oscpar=0);

  void setSeedHisogram(int nselection, TH1D* hseed, int oscpar=0);
  void drawFOM(int nselection, int oscpar);

  private:

  int flgUseEventList;

};

#endif
