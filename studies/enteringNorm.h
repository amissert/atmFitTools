#ifndef NORM_H
#define NORM_H

#define CINTMODE

#include "TH1D.h"
#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "fqProcessedEvent.h"
#include <iostream>
#include "TH2FV.h"


class enteringNorm{
 
 public:

 // const
 enteringNorm(TChain* ch);

 TChain* t2kmc;
 fqProcessedEvent* mcEvent;
 TH1D* hWallNuMu;
 TH1D* hWallNuE;
 TH1D* hWallAll;
 TH1D* hWallPi0;

 void fillWallNuMu();
 void fillWallNuE();
 void fillAll();


};


#ifdef CINTMODE
#include "enteringNorm.cxx"
#endif

#endif
