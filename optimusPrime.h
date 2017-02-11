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
#include "summaryPlots.h"
#include <iostream>
#include <algorithm>

#define PRINTSUMMARY
#define NFOMBINSMAX 100

////////////////////////////////////////
// Class to choose FV optimal point
class optimusPrime{

  public:

  // constructor 
  optimusPrime(TChain* t2kmc,
               int nevents, 
               const char* datadir,
               const char* mapfile);
 
  // vars
  TChain* chmc; //< T2K MC events
  fqProcessedEvent* mcevent; //< individual event values
  int nevents;  //< total number of events to use
  randomList* eventlist; //< random sampling of events if we don't want to use all
  mcLargeArray* fastevents; //< very large array for all T2K FC events
  moreUncertainties* uncertaintyCalculator; //< can calculate addtional uncertainties from entering bg, ect.
  mcmcApply* modifier; //< applies atm fit pars to a given event
  fqcutparams cutPars; //< cut parameters structure defined in eventSelectors.h
  TString mcmcParFileName;
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

  // for calculating figure of merit
//  float Power[NFOMBINSMAX];
//  float Systematic[NFOMBINSMAX];
//  float Nevents[NFOMBINSMAX];

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
  float bestFOMvalue;
  // methods
 
  // masks
  int flgUseMask;
  float maskThresh;
  TH2FV* hMask;
  TH2FV* hMaskOne;
  float smallVariation; 

  
  // use toy-mc to calculate FOM
//  float calcFOMToyMC(float towallmin, float wallmin, int oscpar, int flgselection, int nmcmcpts);


  void fillFVHistoFast(); 
  float calcNuMuFOM(float towallmin, float wallmin, int oscpar);
  float calcFOMSpectrumNuMu(float towallmin, float wallmin, int oscpar, int iplt = 1);
  float calcFOMSpectrumNuE(float towallmin, float wallmin, int oscpar, int iplt = 1);
  float calcFOM(float* pow, float* nev, float* sys, int nbin);
  void calcFOMMap(float towallmax, float wallmax,int oscpar, int npts=15, int flgnumu=1);
  void calcFOMMapE(float towallmax, float wallmax,int oscpar, int npts=15);
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
//  double calcDeltaNuE(float tw1, float w1, float tw2, float w2, int oscpar);

  int isSmallDifference(float tw1, float w1, float tw2, float w2, int oscpar, int flgnumu);

  // superseeds previous methods for applying numu or nue cuts to a modified event
  int applyCutsToModifiedEvent(int iev);
  /////////////////////////////////////////////////////////////////////////////////
 
  void makeAllPlots(float twmax, float wmax, int oscpar, int npts=30,int flgnumu=1);

  void compareCuts(float tw1, float w1, float tw2, float w2, int oscpar, int flgnumu);
  void compareFOM(float tw1, float w1, float tw2, float w2, int oscpar, int flgnumu);
  void showBreakdown();
  void printCutDiff(int flgnumu);
  void printCompare(const char* dir,float tw1, float w1, float tw2, float w2, int oscpar, int flgnumu);

  private:

  int flgUseEventList;

};

#endif
