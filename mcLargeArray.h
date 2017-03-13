#ifndef LARGEARRAY_H
#define LARGEARRAY_H

#define CINTMODE 

#include "TChain.h"
#include "TFile.h"
#include "TMath.h"
#include "TRandom2.h"
#include "fqProcessedEvent.h"

#include <iostream>
#include <algorithm>

#define LARGE 2600000
//#define LARGE 3000000

using namespace std;

class mcLargeArray{

  public:
  // constructor
  mcLargeArray(TChain *ch, int nevents);

  int nsize;

  Short_t vnutype[LARGE];
  UInt_t vnhitac[LARGE];
  float vfqemom[LARGE];
  float vfqmumom[LARGE];
  float vfqrcpar[LARGE];
  Short_t vfqnring[LARGE];
  Short_t vfqnsubev[LARGE];
  Short_t vbin[LARGE];
  Short_t vsample[LARGE];
  Short_t vpassnumu[LARGE];
  Short_t vpassnue1rpi[LARGE];
  Short_t vpassnue[LARGE];
  Short_t vcomponent[LARGE];
  float vfqpi0par[LARGE];
  float vfqpid[LARGE];
  float vfqtowall[LARGE];
  float vfqwall[LARGE];
  Short_t vmode[LARGE];
  float voscpower[LARGE][4];
  float vwallv[LARGE];
  float vtowallv[LARGE];
  float vweight[LARGE];
  float vfqenue[LARGE];
  float vfqenumu[LARGE];
  float vpmomv[LARGE];
  float vattribute[LARGE][5];
  float vfqdir[LARGE][2][3];

  void fillArray(TChain *ch);
  void printEvent(int iev);
//  void fillThinArray(TChain *ch, int nthinning);

};

#ifdef CINTMODE
#include "mcLargeArray.cxx"
#endif

#endif

