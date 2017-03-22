#ifndef MOREUNC_H
#define MOREUNC_H

#include "fqProcessedEvent.h"
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TGraph.h>
#include <TMath.h>
#include <TString.h>
#include <iostream>
#include "TH2FV.h"
#include "TH1D.h"

using namespace std;

///////////////////////////////////////////////////////////////
// fractional uncertaintiy on number of entering backgrounds
// calculated from extrapolating MC
const float enteringNormUncertainty = 0.301;

/////////////////////////////////////////
// uncertainty in the fraction of entering
// events reconstructed inside the ID
// calculated from entering comsics
const float wallNormUncertainty = 0.04;


///////////////////////////////////////
// binning for xsec histos
const int    CCQENBins = 14;
const double CCQEBinning[] = {0.,190.,240.,294.,333.,390.,440.,487.,590.,690.,786.,896.,994.,2000.,3000.};
const double CCQEUncertaintyValues[] = {1.0,0.411,0.216,0.155,0.125,0.105,0.0805,0.055,0.0542,0.0398,0.0344,0.0226,0.0165,0.00903};



///////////////////////////////////////////////////////////////// 
// Class to assign fractional uncertainties to each event.  Includes
// uncertainties from entering backgrounds, FV effects from atmospheric fit.
// This information should be stored int a "data" directory that contines ROOT
// files with the necessary graphs and histograms
class moreUncertainties{

 public:

 // constructor
 moreUncertainties(const char* datadir, const char* mapfile="");
 TString dataDirectory;
 TString mapFileName;

 // graphs and histos of uncertainties
 TGraph* gEnteringWallUnc;
 TGraph* gEnteringNorm;
 TGraph* gEnteringPi0Norm;
 TGraph* gResidualAtt0;
 TGraph* gResidualAtt2;
 TGraph* gResidualAtt3;
 TGraph* gResidualAtt4;
 TH2FV* hfvmap;
 TH2FV* hfvmapccqe;
 TH2FV* hfvmapccnqe;
 TH2FV* hfvmapccwrong;
 TH2FV* hfvmapnc;
 TH1D*  hERecUnc[10]; //< [number of fv bins]
 TH1D*  hCCQEUnc;
 float muMisIDInflation; //< to artificially inflate muon mis-ID background

 // event info
 fqProcessedEvent* mcevent;

 // evaluate errors
 float getTotalUncertainty(float wallv, float wallrc, float towallrc, float erec, int mode, int wronglepton);
 float getXsecUnc(int mode, float enu);
 float getEnteringWallUnc(float wallv, float wallrc);
 float getEnteringWallNormUnc(float wallv);
 float getEnteringNormUnc(float wallv);
 float getFVUncertainty(float towallrc, float wallrc, int mode, int wronglepton);
 float getFVUncEBin(float towallrc, float wallrc, float nuenergy);
 void setEventPointer(fqProcessedEvent* fqevent);

 // initialize
 void init();

 // plottin
 TH2FV* hfv;
 TH1D* hwall;
 TH1D* hwallunc;

};

#ifdef CINTMODE
#include "moreUncertainties.cxx"
#endif

#endif























