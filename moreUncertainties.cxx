#ifndef MOREUNC_C
#define MOREUNC_C

#include "moreUncertainties.h"


moreUncertainties::moreUncertainties(const char* datadir){

 dataDirectory = datadir;

 init();

}

void moreUncertainties::setEventPointer(fqProcessedEvent* fqevent){

 mcevent = fqevent;
 
 return;

}

void moreUncertainties::init(){

  // read in graphs
  TString fname = dataDirectory.Data();
  fname.Append("EnteringUncertainty.root");
  TFile* fileEnteringWallUnc = new TFile(fname.Data());
  gEnteringWallUnc = (TGraph*)fileEnteringWallUnc->Get("entering_bg_uncertainty");

  return;

}

void moreUncertainties::fillFVHisto(TChain* ch){

 TH1D* hwall = new TH1D("hwall","hwall",50,0,200);
 TH1D* hwallunc = new TH1D("hwall","hwall",50,0,200);

 int nev = ch->GetEntries();
 for (int iev=0; iev<nev; iev++){
   ch->GetEntry(iev);
   double wall = mcevent->fqwall;
   double towall = mcevent->fqtowall;
   float totalunc = getTotalUncertainty();
   hwallunc->Fill(wall,mcevent->evtweight*totalunc);
   if (mcevent->wallv<0.) hwall->Fill(wall,mcevent->evtweight);
 }
 
 hwallunc->Divide(hwall);
 hwallunc->Draw("h");

 return;

}

float moreUncertainties::getEnteringNormUnc(){
  if (mcevent->wallv < 0.) return enteringNormUncertainty;
  return 0.;
}

float moreUncertainties::getTotalUncertainty(){

  float totalunc = .0;

  float enteringWallUnc = getEnteringWallUnc();
  totalunc += (enteringWallUnc*enteringWallUnc);

  float normunc = getEnteringNormUnc();
  totalunc += (normunc*normunc);

  return TMath::Sqrt(totalunc);

}


float moreUncertainties::getEnteringWallUnc(){

  // get reconstructed wall
  float  wall = mcevent->fqwall;

  // get true wall
  float wallv = mcevent->wallv;

  if (wallv<=0.) return gEnteringWallUnc->Eval(wall);
  return 0.;

}

#endif




