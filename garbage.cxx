#ifndef HISOTARRAY_H
#define HISOTARRAY_H

#include "modHistoArray.h"



using namespace std;


modHistoArray::modHistoArray(TH1D* hseed,int ninit,const char* tag){

  // get base name tag
  nameTag = hseed->GetName();
  nameTag.Append(tag);

  // setup histogram seed
  TString hname = nameTag.Data();
  hname.Append("_seed");
  hSeed = (TH1D*)hseed->Clone(hname.Data());

  // initialize array
  nHistos = 0;
  currentIndex = 0;
  if (ninit>=0) init(ninit);

}

///////////////////////////////////
//draw all histos
void modHistoArray::drawArray(){
  
  histos[0]->Draw();
  for (int i=0; i<nHistos; i++){
    histos[i]->Draw("same");
  }

  return;
}

///////////////////////////////////
//setup array of histos
void modHistoArray::init(int ninit){

  nHistos = 0;
  for (int i=0; i<ninit; i++){
    TString hname = nameTag.Data();
    hname.Append(Form("_%d",i));
    histos[i] = (TH1D*)hSeed->Clone(hname.Data());
    histos[i]->Reset();
    nHistos++;
  } 
  currentIndex = 0;

  return;
}

///////////////////////////////////////////////////
//copy hisogram contents to a point in the array
void modHistoArray::setHistoContents(TH1D* hadd, int index){

  if (index<0){
    index = currentIndex;
    nHistos++;
  }
 
  int nbins = histos[index]->GetNbinsX();
  for (int ibin=1; ibin<=nbins; ibin++){
    double binc = hadd->GetBinContent(ibin);
    double binerr = hadd->GetBinError(ibin);
    histos[index]->SetBinContent(ibin,binc);
    histos[index]->SetBinError(ibin,binerr);
  }
 
  currentIndex++;
  
  return;
}







#endif

