
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
#include "TGraph.h"


const int NCLASSES = 20;
const int NTOYS    = 500;

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
  float Nevents[NCLASSES][NTOYS];
  float NeventsTotal[NCLASSES][NTOYS];
  float Efficiency[NCLASSES][NTOYS];

  // energy binning
  TH1D* hEvisNuECCQE;
  TH1D* hEvisNuECCOth;
  TH1D* hEvisNuMuCCQE;
  TH1D* hEvisNuMuCCOth;
  // totals
  TH1D* hEvisNuECCQETot;
  TH1D* hEvisNuECCOthTot;
  TH1D* hEvisNuMuCCQETot;
  TH1D* hEvisNuMuCCOthTot;

  TH2D* hCor;
  TH2D* hCov;

  // histogram of all numbers of events
  TH1D* hSlice;
  TGraph* gScat;
//  TH1D* hTmp[NTOYS];

  // initialize histograms
  void initHistos();
  void zeroArrays();
  void resetHistos();

  // draw a particular toy
  void drawSlice(int ntoy);
  void drawSliceEff(int ntoy);
  void drawAll();
  void drawAllEff();

  // for classifying and filling
  int getClassMC(int nutype, int mode, int component, float evis);

  // fill an event in the histograms
//  void addEvent(int nutype, int mode, int component, float evis, float weight); 
  void addEvent(int nclass, float evis, float weight, bool flgtotal); 

  // save histo contents into arrays
  void addToy(int ntoy);

  // calculate correlation and covariance
  void calcCov();

  // calculate correlation and covariance using epsilon
  void calcCovEff();

  // calculate efficiency based on index of total event numbers
  float calcEff(int nclass, int ntoy, int nomindex=0);

  // calculate all effeciencies
  void calcAllEff(int ntoy);;

  // draw scatterplot
  void drawScatter(int iclass, int jclass);

};



#ifdef CINTMODE
#include "SKError.cxx"
#endif










#endif
