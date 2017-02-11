#ifndef OPTIMIZEFV_H
#define OPTIMIZEFV_H

#include "TH2D.h"
#include "TH1D.h"
#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "eventSelector.C"
#include "TH2FV.C"
#include <iostream>
#include "TGraph2D.h"
#include "TGraph.h"


class optimizeFV{
  public:

  // constructor
  optimizeFV();


  // vars
  TChain* chMC;
  fqProcessedEvent* mcEvent

  // methods
  calcFOM();

};





#endif
