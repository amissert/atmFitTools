#ifndef MOREUNC_C
#define MOREUNC_C

#include "moreUncertainties.h"

///////////////////////////////////////////////////////////
// constructor
moreUncertainties::moreUncertainties(const char* datadir){

 dataDirectory = datadir;

 init();

}


///////////////////////////////////////////////////////////////////////
// set the pointer to the event information
void moreUncertainties::setEventPointer(fqProcessedEvent* fqevent){

 mcevent = fqevent;
 
 return;

}

//////////////////////////////////////////////////////
// Get fiducial volume uncertainty
float moreUncertainties::getFVUncertainty(float towallrc, float wallrc, int mode, int wronglepton){

  int ibin = hfvmap->FindBin(towallrc,wallrc);
  if (ibin<=0) return 0;
  float sys=0.;
  if (wronglepton){
    sys = hfvmapccwrong->GetBinContent(ibin);
  }
  else if (mode*mode==1){
    sys = hfvmapccqe->GetBinContent(ibin); 
  }
  else if (TMath::Abs(mode)<30){
    sys = hfvmapccnqe->GetBinContent(ibin);
  }
  else{
    sys = hfvmapnc->GetBinContent(ibin);
  }
  

//  if (ibin<=0) return 0;
//  else{
//    float syst = hfvmap->GetBinContent(ibin);
//
    return sys;
//  }
}


///////////////////////////////////////////////////////
// get the fractional uncertainty in specfice Erec bin
// in a specific FV area
float moreUncertainties::getFVUncEBin(float towallrc, float wallrc, float erec){

  // get FV bin
  int ibin = hfvmap->FindBin(towallrc,wallrc);

  // if outside FV, do nothing
  if (ibin<=0) return 0;
  if (ibin==6) return 0;

  // otherwise, get the uncertainty in this bin
  int ebin = hERecUnc[ibin-1]->FindBin(erec);
  return (float) hERecUnc[ibin-1]->GetBinContent(ebin);

}


///////////////////////////////////////////////////////
// Read in some important graphs and histograms from the 
// data directory
void moreUncertainties::init(){

  // read in graphs
  TString fname = dataDirectory.Data();
  fname.Append("EnteringUncertainty.root");
  TFile* fileEnteringWallUnc = new TFile(fname.Data());
  gEnteringWallUnc = (TGraph*)fileEnteringWallUnc->Get("entering_bg_uncertainty");

  // read in FV uncertainty maps
  TString fmapname = dataDirectory.Data();
  fmapname.Append("FVUncMap.root");
  TFile *filefvmap = new TFile(fmapname.Data());
  hfvmap = (TH2FV*)filefvmap->Get("FVUncMapNC");
  hfvmapccqe = (TH2FV*)filefvmap->Get("FVUncMapCCQE");
  hfvmapccnqe = (TH2FV*)filefvmap->Get("FVUncMapCCnQE");
  hfvmapccwrong = (TH2FV*)filefvmap->Get("FVUncMapCCWrong");
  hfvmapnc = (TH2FV*)filefvmap->Get("FVUncMapNC");
  
  // read in uncertainties for Erec bins
  for (int ifvbin=0; ifvbin<hfvmap->GetNumberOfBins(); ifvbin++){
    hERecUnc[ifvbin] = (TH1D*)filefvmap->Get(Form("Bin_Uncertainty_FVBin%d",ifvbin));
  }

  //
  return;

}

/////////////////////////////////////////////////////////
// Fills a histogram of uncertaintes, can be useful for
// debugging
/*
void moreUncertainties::fillFVHisto(TChain* ch){

  hwall = new TH1D("hwall","hwall",50,0,200);
  hwallunc = new TH1D("hwall","hwall",50,0,200);

 int nev = ch->GetEntries();
 for (int iev=0; iev<nev; iev++){
   ch->GetEntry(iev);
   double wall = mcevent->fqwall;
   double towall = mcevent->fqtowall;
   float totalunc = getTotalUncertainty(mcevent->wallv,wall,towall,mcevent->fq1renu[1]);
   hwallunc->Fill(wall,mcevent->evtweight*totalunc);
   if (mcevent->wallv<0.) hwall->Fill(wall,mcevent->evtweight);
 }
 
 hwallunc->Divide(hwall);
 hwallunc->Draw("h");

 return;

}
*/

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


////////////////////////////////////////////////////////////////
// get the x-section uncertainty
float moreUncertainties::getXsecUnc(int mode){

  int absmode = TMath::Abs(mode);

  if (absmode==1){
    return 0.07;
  }
  else if (absmode<30){
    if (absmode==16){
      return 1.0;
    }
    else{
      return 0.2;
    }
  }
  else if (absmode==36){
    return 1.0;
  }
  return 0.3;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Return the total fractional uncertainty for this event
float moreUncertainties::getTotalUncertainty(float wallv, float wallrc, float towallrc, float erec, int mode, int wronglepton ){

  
  float totalunc = .0;

  float enteringWallUnc = getEnteringWallUnc(wallv,wallrc);
  totalunc += (enteringWallUnc*enteringWallUnc);

  float normunc = getEnteringNormUnc(wallv);
  totalunc += (normunc*normunc);

  float wallnormunc = getEnteringWallNormUnc(wallv);
  totalunc += (wallnormunc*wallnormunc);

  // for FV
//  float fvuncebin = getFVUncEBin(towallrc,wallrc,erec);
//  totalunc += (fvuncebin*fvuncebin);
 
  // for x section
  float xunc = getXsecUnc(mode);
  totalunc += xunc;

//     for FV
  float fvunc = getFVUncertainty(towallrc,wallrc,mode,wronglepton);
  totalunc += (fvunc*fvunc);

  float syst = TMath::Sqrt(totalunc);
  
//  cout<<"totalunc: "<<syst<<endl;
  return syst;
//  return 0.1;
}

///////////////////////////////////////////////////
// fractional uncertainty from wall shape
float moreUncertainties::getEnteringWallUnc(float wallv, float wallrc){

  if (wallv<=0.) return gEnteringWallUnc->Eval(wallrc);

  return 0.;

}

#endif




