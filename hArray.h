
#ifndef HARRAY_H
#define HARRAY_H

#include "TH1D.h"

#include "TH2D.h"
#include "TString.h"
#include "stats.C"


#define HARRAYMAX 500

class hArray{

  public:

  hArray(const char* name, TH1D* hseed, int nh);

  TH1D* histos[HARRAYMAX];
  TString nameTag;
  TH1D* hSeed;
  TH1D* hSummary;
  int nHistos;

  void init();
  void calcSummary();
  void setHistogram(int ih, TH1D* hh);
  void drawArray();
  void delArray();
//  void cm
};

#include "hArray.cxx"

#endif
