#ifndef SKERROR_C
#define SKERROR_C

#include "SKError.h"



///////////////////////////////////////////////////////////////
void SKError::addToy(int ntoy){

 // save histos in arrays
 int index = 0;
 //
 for (int ibin=1; ibin<=hEvisNuECCQE->GetNbinsX(); ibin++){
   Nevents[index][ntoy] = hEvisNuECCQE->GetBinContent(ibin);
   index++;
 }
 //
 for (int ibin=1; ibin<=hEvisNuECCOth->GetNbinsX(); ibin++){
   Nevents[index][ntoy] = hEvisNuECCOth->GetBinContent(ibin);
   index++;
 }
 //
 for (int ibin=1; ibin<=hEvisNuMuCCQE->GetNbinsX(); ibin++){
   Nevents[index][ntoy] = hEvisNuMuCCQE->GetBinContent(ibin);
   index++;
 }
 //
 for (int ibin=1; ibin<=hEvisNuMuCCOth->GetNbinsX(); ibin++){
   Nevents[index][ntoy] = hEvisNuMuCCOth->GetBinContent(ibin);
   index++;
 }

 // clear histos for next toy
 resetHistos();

 // 
 return;

}


///////////////////////////////////////////////////////////////
void SKError::resetHistos(){

  //
  hEvisNuECCQE->Reset();
  hEvisNuMuCCQE->Reset();
  hEvisNuECCOth->Reset();
  hEvisNuMuCCOth->Reset();

  //
  return;
}


///////////////////////////////////////////////////////////////
void SKError::addEvent(int nclass, float evis, float weight){

  //
  if (nclass==1){
    hEvisNuECCQE->Fill(evis,weight);
  }
  if (nclass==2){
    hEvisNuMuCCQE->Fill(evis,weight);
  }
  if (nclass==3){
    hEvisNuECCOth->Fill(evis,weight);
  }
  if (nclass==4){
    hEvisNuMuCCOth->Fill(evis,weight);
  }

  //
  return;
}


///////////////////////////////////////////////////////////////
int SKError::getClass(int nutype, int mode, int component, float evis){
  
  // class code: 
  //  1 -> CCQE Nu E
  //  2 -> CCQE Nu Mu
  //  3 -> CCOth Nu E
  //  4 -> CCOth Nu Mu
 
  if (component==0){
    return 1;
  }
  if (component==2){
    return 3;
  }
  if (component==1){
    return 2;
  }
  if (component==3){
    return 4;
  }

  return 0;
}



///////////////////////////////////////////////////////////////
void SKError::drawSlice(int ntoy){

  // reset and fill 
  hSlice->Reset();
  for (int iclass=0; iclass<Nclass; iclass++){
    hSlice->SetBinContent(iclass+1,Nevents[iclass][ntoy]);
    hSlice->SetBinError(iclass+1,0.);
  }

  hSlice->Draw();
  
  //
  return; 
}



///////////////////////////////////////////////////////////////
void SKError::zeroArrays(){

  // set initial values
  for (int i=0; i<Nclass; i++){
    Nevents_nominal[i] = 0.;
    for (int j=0; j<Nclass; j++){
      Nevents[i][j] = 0.;
    }
  }

  return;
}



///////////////////////////////////////////////////////////////
void SKError::initHistos(){

  // histo binning 
  int    NbinsNuECCQE = 6;
  double BinsNuECCQE[] = {100.,300.,700.,1250.,2000.,5000.,30000.};
  int    NbinsNuECCOth = 3;
  double BinsNuECCOth[] = {100.,1250.,5000.,30000.};
  int    NbinsNuMuCCQE = 7;
  double BinsNuMuCCQE[] = {0.,100.,300.,700.,1250.,2000.,5000.,30000.};
  int    NbinsNuMuCCOth = 3;
  double BinsNuMuCCOth[] = {100.,1250.,5000.,30000.};

  // setup histos
  hEvisNuECCQE = new TH1D("hnueccqe","hnueccqe",NbinsNuECCQE,BinsNuECCQE); 
  hEvisNuECCOth = new TH1D("hnueccqe","hnueccqe",NbinsNuECCOth,BinsNuECCOth); 
  hEvisNuMuCCQE = new TH1D("hnumuccqe","hnumuccqe",NbinsNuMuCCQE,BinsNuMuCCQE); 
  hEvisNuMuCCOth = new TH1D("hnumuccoth","hnumuccoth",NbinsNuMuCCOth,BinsNuMuCCOth); 

  // count all bins
  Nclass = NbinsNuECCQE + NbinsNuECCOth + NbinsNuMuCCQE + NbinsNuMuCCOth;

  // setup slice
  hSlice = new TH1D("hslice","hslice",Nclass,0,Nclass);

  //
  return;
}



///////////////////////////////////////////////////////////////
SKError::SKError(int ntoys){

  // setup the histos
  Ntoys = ntoys;
  initHistos();

}





#endif
