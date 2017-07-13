#ifndef FITPLOTS_H
#define FITPLOTS_H

#include "histoCompare.h"
#include "atmFitPars.h"
#include "hArray.h"
#include "TRandom2.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "mcmcReader.h"

using namespace std;

//////////////////////////////////////////////
// loop over ntestpts
// for each pt:
//   set atmFitPars from that point
//   get modified histogram
//   add mod histogram to arrays
//
class fitPlots{

  public:
 
  // constructor
  fitPlots(histoCompare* hc, TTree* partree);

  // vars
  TTree* parTree;
  histoCompare* hCompare;
  mcmcReader* mcmcPars;
  hArray* hAtt[5];
  int nAtt;
  int nPoints;
  TCanvas* cc;
  TH1D* hTmp;

  // methods
  void applyPars();
  void fillArrays(int ibin, int isamp);
  void drawFitSummary(int isamp, int ibin); //< draw all attributes
  void drawFitSummaryFV(int isamp, int iatt); //< draw all bins
  void drawFitThrows(int isamp, int ibin);
  void printFitSummary(const char* dir, const char* name);
  void printFitSummaryFV(const char* dir, const char* name);
  void initArrays();

  TString getAxisTitle(int iatt);

};


#endif
