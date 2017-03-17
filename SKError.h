
#ifndef SKERROR_H
#define SKERROR_H



#include <iostream>

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
//
class SKError{

  public:

  // constructor
  SKError(int ntoy = 100);
  int Nclass;
  int Ntoys;
 
  // arrays for numbers of events
  double Nevents[NCLASSES][NTOYS];
  double NeventsTotal[NCLASSES][NTOYS];
//  double Efficiency[NCLASSES][NTOYS];
  double DelEfficiency[NCLASSES][NTOYS];
  double DelEffShiftError[NCLASSES];
  double absDelEffShiftError[NCLASSES];
  double DelEffFitError[NCLASSES];


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


  // for classifying and filling
  int getClassMC(int nutype, int mode, int component, double evis, int nsubev, double towall, double wall);

  // fill an event in the histograms
  void addEvent(int nclass, double evis, double weight, bool flgtotal); 

  // save histo contents into arrays
  void addToy(int ntoy);

  // calculate correlation and covariance
  void calcCov();

  // calculate correlation and covariance using epsilon
  void calcCovEff();

  // calculate correlation and covariance using epsilon
  void calcCovDelEff();

  // calculate efficiency based on index of total event numbers
  double calcEff(int nclass, int ntoy);

  // calculate efficiency based on index of total event numbers
  double calcDelEff(int nclass, int ntoy);

  // calculate all effeciencies
  void calcAllDelEff(int ntoy);

  // draw scatterplot
  void drawScatter(int iclass, int jclass);

  //
  double calcShiftError(int iclass);

  //
  double calcFitError(int iclass);

  //
  void calcErrors();

  // 
  void printErrors();

  void saveErrors(const char* filename);
  

  void printEffDist(const char* plotdir);

};



#ifdef CINTMODE
#include "SKError.cxx"
#endif










#endif
