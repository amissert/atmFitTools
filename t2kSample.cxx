#ifndef T2KSAMPLE_CXX
#define T2KSAMPLE_CXX

#include "t2kSample.h"



///////////////////////////////////////////////////
// calculate errors
void t2kSample::calcUncertainties(){

  return;
}


///////////////////////////////////////////////////
// read all histos
void t2kSample::readHistos(){

  return;
}


///////////////////////////////////////////////////
// save all histos
void t2kSample::saveHistos(const char* dir){

//  TString filename = dir;
//  filename.Append("T2K_Selection");
//  filename.Append(nameTag.Data());
//  filename.Append(".root");

//  TFile* fout = new TFile(filename.Data(),"RECREATE");
//  fout->cd();
  
  
  return;

}


///////////////////////////////////////////////////
// update array contents after finishing a toy run
void t2kSample::finishToyRun(){
  
  // update arrays with bin contents
  int iHistoBin = 0;

  if (flgE){
    for (int i=0; i<=hEnuElectron->GetNbinsX(); i++){
      enuBinValues[toyCounter][iHistoBin] = (float)hEnuElectron->GetBinContent(i);
      iHistoBin++;
    }
  }

  if (flgMu){
    for (int i=0; i<=hEnuMuon->GetNbinsX(); i++){
      enuBinValues[toyCounter][iHistoBin] = (float)hEnuMuon->GetBinContent(i);
      iHistoBin++;
    }
  }

  resetAll();
  toyCounter++;

  //
  return;
}


//////////////////////////////////////////
// reset values between toys
void t2kSample::resetAll(){
  hEnuElectron->Reset();
  hEnuElectron->Reset();
  
  //
  return;
}


//////////////////////////////////////////
// initialize histograms
void t2kSample::initHistos(){

  TString hname = "ElectronSample";
  hname.Append(nameTag.Data());
  hEnuElectron = new TH1D(hname.Data(),hname.Data(),NEnuBinsElectrons,EnuBinningElectrons);
  hEnuMuon = new TH1D(hname.Data(),hname.Data(),NEnuBinsMuons,EnuBinningMuons);

  nTotBins = 0;
  if (flgE) nTotBins += hEnuElectron->GetNbinsX();
  if (flgMu) nTotBins += hEnuMuon->GetNbinsY();

 
  nTotBins = hEnuElectron->GetNbinsX() + hEnuMuon->GetNbinsX(); 

  //
  return;
}


//////////////////////////////////////////
// constructor
t2kSample::t2kSample(const char* nametag, int flgelectron, int flgmuon){

  nameTag = nametag;
  flgE = flgelectron;
  flgMu = flgmuon;
  toyCounter = 0;
  initHistos();

}



#endif
