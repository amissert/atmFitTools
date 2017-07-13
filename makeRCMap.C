#ifndef  makeRCMap_INC
#define  makeRCMap_INC

#include <iostream>

#include "TString.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"

#include "fqProcessedEvent.cxx"

//#define CINTMODE

/////////////////////////////////////////////////////////////
class makeRCMap{
  public:

  makeRCMap(const char* files);

  TChain* ch;

  void makeMap();
  TH2D* map;

  void init();
};



////////////////////////////////////////////////////////////
void makeRCMap::makeMap(){
  int nev = ch->GetEntries();
  fqProcessedEvent* fqevent = new fqProcessedEvent(ch);
  for (int iev=0; iev<nev; iev++){
    ch->GetEntry(iev);
    map->Fill(fqevent->fqmrnring[0],fqevent->attribute[3],fqevent->evtweight);
  }
  map->Draw("colz");
  return;
}



////////////////////////////////////////////////////////////
void makeRCMap::init(){

  int nbiny = 500;
  int nbinx = 7;
  double xmin = 0;
  double xmax = 7;
  double ymin = -2000;
  double ymax = 2000;
  map = new TH2D("hmap","hmap",nbinx,xmin,xmax,nbiny,ymin,ymax);
  return;
}



////////////////////////////////////////////////////////////
makeRCMap::makeRCMap(const char* files){
  ch = new TChain("h1");
  ch->Add(files);
  init();
}


#endif   /* ----- #ifndef makeRCMap_INC  ----- */




