#ifndef  RCMap_H
#define  RCMap_H

#include <iostream>

#include "TMath.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"


///////////////////////////////////////////////////
class RCMap{
  public:
  RCMap(const char* fname);

  TH2D* hmap;
  TH2D* hweight;

  void init();
  double getRowSum(TH2D* h, int irow);
  void scaleRow(TH2D* h, int irow, double scale);
  void getConditionalX();
  double get1RWeight(double rcpar);

};

#ifdef CINTMODE
#include "RCMap.C"
#endif

#endif
