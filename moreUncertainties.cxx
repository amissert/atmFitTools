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

  hwall = new TH1D("hwall","hwall",50,0,200);
  hwallunc = new TH1D("hwall","hwall",50,0,200);

 int nev = ch->GetEntries();
 for (int iev=0; iev<nev; iev++){
   ch->GetEntry(iev);
   double wall = mcevent->fqwall;
   double towall = mcevent->fqtowall;
   float totalunc = getTotalUncertainty(mcevent->wallv,wall);
   hwallunc->Fill(wall,mcevent->evtweight*totalunc);
   if (mcevent->wallv<0.) hwall->Fill(wall,mcevent->evtweight);
 }
 
 hwallunc->Divide(hwall);
 hwallunc->Draw("h");

 return;

}


//////////////////////////////////////////////////////
// fractional uncertainty from not simulating all of dead region
float moreUncertainties::getEnteringWallNormUnc(float wallv){
  if (wallv < 0.) return wallNormUncertainty;
  return 0.;
}


//////////////////////////////////////////////////////
// fractional uncertainty from not simulating all of dead region
float moreUncertainties::getEnteringNormUnc(float wallv){
  if (wallv < 0.) return enteringNormUncertainty;
  return 0.;
}

float moreUncertainties::getTotalUncertainty(float wallv, float wallrc){

  float totalunc = .0;

  float enteringWallUnc = getEnteringWallUnc(wallv,wallrc);
  totalunc += (enteringWallUnc*enteringWallUnc);

  float normunc = getEnteringNormUnc(wallv);
  totalunc += (normunc*normunc);

  float wallnormunc = getEnteringWallNormUnc(wallv);
  totalunc += (wallnormunc*wallnormunc);

  return TMath::Sqrt(totalunc);

}

///////////////////////////////////////////////////
// fractional uncertainty from wall shape
float moreUncertainties::getEnteringWallUnc(float wallv, float wallrc){

  if (wallv<=0.) return gEnteringWallUnc->Eval(wallrc);

  return 0.;

}

#endif




