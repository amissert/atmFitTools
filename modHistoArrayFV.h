#ifndef HISTOARRAYFV_H
#define HISTOARRAYFV_H

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


using namespace std;


#define ARRAYMAX 500
#define FVBINSMAX 10

/////////////////////////////////////////////////////////////
// A class to manage arrays histograms, where each histogram
// is filled using some systematic parameter throw. 
// Construct by seeding with a 1D histogram and a TH2FV.  
// An array of 1D histogram will be created, with an idex corresponding
// to the TH2FV bin and an index corresponding to the systematics throw.
class modHistoArrayFV{

  public:

  // constructor
  modHistoArrayFV(TH1D* hseed, TH2FV* hfv, int ninit=-1);
  modHistoArrayFV(const char* filename, int nfvbins, int nthrows);

  // initialize 1D histograms
  void init(int ninit=-1);

  // internal vars
  TFile* fout;
  TH1D* histos[ARRAYMAX][FVBINSMAX];
  int nHistos; //< total number
  int nPoints;
  TH1D* hSeed;
  TH2FV* hFV[ARRAYMAX];
  TH2FV* hFVCCQE[ARRAYMAX];
  TH2FV* hFVCCnQE[ARRAYMAX];
  TH2FV* hFVCCWrong[ARRAYMAX];
  TH2FV* hFVNC[ARRAYMAX];
  TString nameTag;
  int currentIndex;

  // summary histograms
  // bin uncertainties
  TH1D*  binUnc[FVBINSMAX]; 
  TLine* nominalLine[FVBINSMAX];
  //
  TH2FV* FVUncMap; //< total uncertainty in # of events
  TH2FV* FVUncMapCCQE; //< total uncertainty in # of events
  TH2FV* FVUncMapCCnQE; //< total uncertainty in # of events
  TH2FV* FVUncMapCCWrong; //< total uncertainty in # of events
  TH2FV* FVUncMapNC; //< total uncertainty in # of events
  //
  TH2FV* FVShiftMap;
  TH2FV* FVShiftMapCCQE;
  TH2FV* FVShiftMapCCnQE;
  TH2FV* FVShiftMapCCWrong;
  TH2FV* FVShiftMapNC;
  //
  TH2FV* FVFitMap;
  TH2FV* FVFitMapCCQE;
  TH2FV* FVFitMapCCnQE;
  TH2FV* FVFitMapCCWrong;
  TH2FV* FVFitMapNC;

  TH1D*  hNevents[FVBINSMAX];
  TH1D*  hNeventsCCQE[FVBINSMAX];
  TH1D*  hNeventsCCnQE[FVBINSMAX];
  TH1D*  hNeventsCCWrong[FVBINSMAX];
  TH1D*  hNeventsNC[FVBINSMAX];

  // methods
  void readFromFile(int nfvbins, int nthrows); //< read in previously filld histograms
  TH1D* getHistogram(int index, int fvbin); //< return a particluar 1D histogram from an FV region
  void setHistoContents(TH1D* hadd, int index, int fvbin); 
  void drawArray(int fvbin);
  void drawNev(int fvbin);
  void calcSummary();
  void saveSummary(const char* fname);
  void saveClose();
//  void makeUncMap(const char* outfile);
  void printUncMap(const char* plotdir);
//  void printUncSummary(const char* plotdir, int type);
  private:
  void fillFitAndShiftError();


};


#ifdef CINTMODE
#include "modHistoArrayFV.cxx"
#endif

#endif
