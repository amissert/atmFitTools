#ifndef MOREUNC_C
#define MOREUNC_C

#include "moreUncertainties.h"

///////////////////////////////////////////////////////////
// constructor
moreUncertainties::moreUncertainties(const char* datadir, const char* mapfile){

 dataDirectory = datadir;
 mapFileName = mapfile;

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

  // find which FV bin we are in
  int ibin = hfvmap->FindBin(towallrc,wallrc);
  if (ibin<=0) return 0;
  float sys=0.;
  
  // is event from NC?
  if (TMath::Abs(mode)>=30){
    sys = hfvmapnc->GetBinContent(ibin);
    return sys;
  }
  // is event Mis-ID?
  else if (wronglepton){
    sys = hfvmapccwrong->GetBinContent(ibin) + muMisIDInflation;
    return sys;
  }
  // now for ccqe
  else if (TMath::Abs(mode)<=10){
    sys = hfvmapccqe->GetBinContent(ibin); 
    return sys;
  }
  // for cc nqe
  else if (TMath::Abs(mode)<30){
    sys = hfvmapccnqe->GetBinContent(ibin);
    return sys;
  }
  

  //
  return sys;
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

  // default to nu-mu map
  if (!mapFileName.CompareTo("")){
    mapFileName = "FVUncMapNuMu.root";
  }
   
  // open map file and read in histograms 
  fmapname.Append(mapFileName.Data());
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

  // for xsec unc
  hCCQEUnc = new TH1D("hccqeunc","hccqeunc",CCQENBins,CCQEBinning);
  for (int i=1; i<= CCQENBins; i++){
    hCCQEUnc->SetBinContent(i,CCQEUncertaintyValues[i-1]);
  }

  // set artificial misID inflation to zero
  muMisIDInflation = 0.;

  //
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



////////////////////////////////////////////////////////////////
// get the x-section uncertainty
float moreUncertainties::getXsecUnc(int mode, float enu){

  int absmode = TMath::Abs(mode);

  // CCQE
  if (absmode<=10){
    int ibin = hCCQEUnc->FindBin(enu);
    if (ibin>0){
       return hCCQEUnc->GetBinContent(ibin);
    }
    else{
      return 0.;
    }
  }
  
  // CCnQE
  else if (absmode<30){
    return 0.2;
  }
 
  // NC Coh
  else if (absmode==36){
    return 1.0;
  }

  // NC
  return 0.2;
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Return the total fractional uncertainty for this event
float moreUncertainties::getTotalUncertainty(float wallv, float wallrc, float towallrc, float erec, int mode, int wronglepton ){

  
  float totalunc = .0;

  // entering from wall shape uncertainty
  float enteringWallUnc = getEnteringWallUnc(wallv,wallrc);
  totalunc += (enteringWallUnc*enteringWallUnc);

  // normalization uncertainty for entering bg
  float normunc = getEnteringNormUnc(wallv);
  totalunc += (normunc*normunc);

  // additional normalization uncertainty from fraction of entering events
  float wallnormunc = getEnteringWallNormUnc(wallv);
  totalunc += (wallnormunc*wallnormunc);
 
  // for x section
  float xunc = getXsecUnc(mode, erec);
  totalunc += xunc*xunc;

  //  for FV
  float fvunc = getFVUncertainty(towallrc,wallrc,mode,wronglepton);
  totalunc += (fvunc*fvunc);

  float syst = TMath::Sqrt(totalunc);
  
  return syst;
}



///////////////////////////////////////////////////
// fractional uncertainty from wall shape
float moreUncertainties::getEnteringWallUnc(float wallv, float wallrc){

  if (wallv<=0.) return gEnteringWallUnc->Eval(wallrc);

  return 0.;

}



#endif




