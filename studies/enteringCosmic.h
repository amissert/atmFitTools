#ifndef COSMIC_H
#define COSMIC_H

#define CINTMODE

#include "TH1D.h"
#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TGraph.h"
#include <iostream>
#include "fqProcessedEvent.h"

class enteringCosmic{

 public:

 enteringCosmic(TChain* chdata, TChain* chmc);

 TChain *chMC;
 TChain *chData;
 fqProcessedEvent* mcevent;
 fqProcessedEvent* datevent;

 double calcNorm;
 double Norm;

 void fillFVHistos();
 void makeGraphs();
 void init();
 TH1D* hwall[2];
 TH1D* hwalle[2];
 TH1D* hwallraw[2];
 TH1D* hratio;
 TH1D* hratioe;
 TH1D* hratiosmooth;
 TH1D* hratioesmooth;
 TH1D* hscale;
 TGraph* gUncMu;
 TGraph* gUncE;
 
 
};


#endif
