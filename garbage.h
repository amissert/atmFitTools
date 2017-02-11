#ifndef HISTOARRAY_H
#define HISTOARRAY_H

#include "TH1D.h"

#include "TH2D.h"
#include "TString.h"


using namespace std;


#define HISTARRAYMAX 50

///////////////////////////////////////////////////
// A class to manage arrays of similar histograms
class modHistoArray{

  public:

  // constructor
  modHistoArray(TH1D* hseed, int ninit=-1,const char* tag="");

  // initialized
  void init(int ninit=-1);

  // internal vars
  TH1D* histos[HISTARRAYMAX];
  int nHistos;
  TH1D* hSeed;
  TString nameTag;
  int currentIndex;

  // methods
  void setHistoContents(TH1D* hadd,int index=-1);
  void drawArray();

};


#ifdef CINTMODE
#include "modHistoArray.cxx"
#endif

#endif
