
#ifndef SKERROR_H
#define SKERROR_H



#include <iostream>
#include <vector>
#include <cstdio>

#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TH2D.h"
#include "TChain.h"
#include "TMath.h"
#include "TString.h"
#include "stats.C"
#include "TLine.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "mcLargeArray.h"
#include "TLatex.h"

const int NCLASSES = 50;
const int NTOYS    = 2000;
const int NLINES= 3;

using namespace std;


///////////////////////////////////////////////////////////////////
class SKError{

  public:

  // constructor
  SKError(int ntoy = 100);
  int Nclass;
  int Ntoys;
  int Nmarginal;
  int NEvisHistos;
  bool flgNC;
 
  // arrays for numbers of events
  double Nevents[NCLASSES][NTOYS];
  double NeventsNominal[NCLASSES];
  double NeventsNominalCore[NCLASSES];
  double NeventsTotal[NCLASSES][NTOYS];
  double DelEfficiency[NCLASSES][NTOYS];
  double DelEffShiftError[NCLASSES];
  double absDelEffShiftError[NCLASSES];
  double DelEffFitError[NCLASSES];

  // compare to TN186 (Table 8)
  double tn186ShiftError[NCLASSES];
  double tn186FitError[NCLASSES];
  double tn186TotError[NCLASSES];
  TH1D* hErrorTN186CCQE[2];
  TH1D* hErrorTN186CCOth[2];

  // evis binning
  TH1D* hEvisNuECCQE;
  TH1D* hEvisNuECCOth;
  TH1D* hEvisNuMuCCQE;
  TH1D* hEvisNuMuCCOth;
  TH1D* hEvisNuENCPi0;
  TH1D* hEvisNuENCPi;
  TH1D* hEvisNuMuNCPi0;
  TH1D* hEvisNuMuNCPi;
  TH1D* hEvisCore[NCLASSES];

  // totals
  TH1D* hEvisNuECCQETot;
  TH1D* hEvisNuECCOthTot;
  TH1D* hEvisNuMuCCQETot;
  TH1D* hEvisNuMuCCOthTot;
  TH1D* hEvisNuENCPi0Tot;
  TH1D* hEvisNuENCPiTot;
  TH1D* hEvisNuMuNCPi0Tot;
  TH1D* hEvisNuMuNCPiTot;
  TH1D* hEvisBase[NCLASSES];

  // for checking distributions
  TH1D* hVariable1;
  TH1D* hVariable2;

  // names and energy ranges
  TString className[NCLASSES];
  TString classRange[NCLASSES];

  // for drawing lines
  TLine* lineHorz[NLINES];
  TLine* lineVert[NLINES];
  vector<TLatex*> vLabels;
  TLatex* labelVert[NCLASSES];
  TLatex* labelHorz[NCLASSES];
  TLatex* sectorLabelHorz[4];
  TLatex* sectorLabelVert[4];
  TLatex* nuLabelVert[4];
  TLatex* nuLabelHorz[4];
  int    lineVal[NLINES];  

  TH2D* hCor;
  TH2D* hCov;
  TVectorD* vShiftErrors;
  TH1D* hErrors;
  TH1D* hErrorsCCQE[2];
  TH1D* hErrorsCCOth[2];

  // histogram of all numbers of events
  TH1D* hSlice;
  TGraph* gScat;
  TH1D* hdist;
  TLine* distMean;

  // flags etc.
  int effDefinition;

  // initialize histograms
  void initHistos(int ibinning=0);
  void zeroArrays();
  void resetHistos();

  // draw a particular toy
  void drawSlice(int ntoy);
  void drawSliceTot(int ntoy);
  void drawSliceEff(int ntoy);
  void drawAll();
  void drawAllEff();
  void drawDist(int nclass);
  void drawEffDist(int nclass);
  void drawCor();
  void drawCov();
  void drawErrors1D(bool flgDrawTN186=false);
  void calcErrors1D();


  // for classifying and filling
  int getClassMC(int nutype, int mode, int component,
                 double evis, int nsubev, double towall, double wall);

  int getClassEvis(int nclass, int nevisbin);

  // fill an event in the histograms
  int addEvent(int nclass, double evis, double weight, bool flgtotal); 

  // fill an event in the histograms
  int addEventBase(int nclass, double evis, double weight); 

  // fill an event in the histograms
  int addEventCore(int nclass, double evis, double weight, int coretype=1); 

  // save histo contents into arrays
  void addToy(int ntoy);
  
  // save histo contents into arrays
  void addNominal();

  // calculate correlation and covariance
  void calcCov(int vartype=0);

  // calculate efficiency based on index of total event numbers
  double calcEff(int nclass, int ntoy);

  // calculate efficiency based on index of total event numbers
  double calcDelEff(int nclass, int ntoy);

  // calculate efficiency based on index of total event numbers
  double calcMargEff(int nclass);

  // calculate efficiencies for all classes
  void calcAllEff(int ntoy);

  // fill Nevent arrays from histograms
  void fillNeventArrays(int ntoy);

  // draw scatterplot
  void drawScatter(int iclass, int jclass);

  // calc the shift error from nominal
  double calcShiftError(int iclass);

  // calc error from variance
  double calcFitError(int iclass);

  // calc all of the fit and shift errors
  void calcErrors();

  // calcuate marginalized efficiency and set as ntoy-th 
  // array value of DelEfficiency
  void marginalize(int ntoy);

  // print out error values
  void printErrors();

  void saveErrors(const char* filename);

  void printEachEffDist(const char* plotdir);

  void printAllEffDist(const char* plotdir, int nclass, const char* tag);

  void makeBinLabels();

  void drawBinLabels();

  void drawVertLines();

  void drawHorizLines();

  vector<TLatex*> getBinLabels(TH1D* hh);

  double getMaxError(TH1D* hh);

  void initTN186Errors();

  void initClassTags();

  void drawVariable(mcLargeArray* events, const char* varname, int iclass, int flgcore=0);

  TH1D* makeHisto(const char* varname, const char* hname);

  int isCore(int nclass, float fqpidpar, float fqpi0par, float fqpippar,float fqrcpar);

  private:

  int evisClass[NCLASSES][NCLASSES];
  
};



#ifdef CINTMODE
#include "SKError.cxx"
#endif










#endif
