
#ifndef SKERROR_H
#define SKERROR_H



#include <iostream>
#include <vector>
#include <cstdio>

#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TChain.h"
#include "TMath.h"
#include "TString.h"
#include "stats.C"
#include "TLine.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TGraph.h"
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
 
  // arrays for numbers of events
  double Nevents[NCLASSES][NTOYS];
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

  // totals
  TH1D* hEvisNuECCQETot;
  TH1D* hEvisNuECCOthTot;
  TH1D* hEvisNuMuCCQETot;
  TH1D* hEvisNuMuCCOthTot;

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
  TH1D* hDiagonalErrors;
  TH1D* hDiagonalErrorsCCQE[2];
  TH1D* hDiagonalErrorsCCOth[2];

  // histogram of all numbers of events
  TH1D* hSlice;
  TGraph* gScat;
  TH1D* hdist;
  TLine* distMean;
  TLine* zeroValue;


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
  void drawDiagonals();
  void calcDiagonals();


  // for classifying and filling
  int getClassMC(int nutype, int mode, int component,
                 double evis, int nsubev, double towall, double wall);

  // fill an event in the histograms
  void addEvent(int nclass, double evis, double weight, bool flgtotal); 

  // save histo contents into arrays
  void addToy(int ntoy);

  // calculate correlation and covariance
  void calcCov(int vartype=0);

  // calculate efficiency based on index of total event numbers
  double calcEff(int nclass, int ntoy);

  // calculate efficiency based on index of total event numbers
  double calcDelEff(int nclass, int ntoy);

  // calculate all effeciencies
  void calcAllDelEff(int ntoy, int effdef=0);

  // draw scatterplot
  void drawScatter(int iclass, int jclass);

  double calcShiftError(int iclass);

  double calcFitError(int iclass);

  void calcErrors();

  void printErrors();

  void saveErrors(const char* filename);

  void printEffDist(const char* plotdir);

  void makeBinLabels();

  void drawBinLabels();

  void drawVertLines();

  void drawHorizLines();

  vector<TLatex*> getBinLabels(TH1D* hh);

  double getMaxError(TH1D* hh);

  void initTN186Errors();


};



#ifdef CINTMODE
#include "SKError.cxx"
#endif










#endif
