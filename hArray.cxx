
#ifndef HARRAY_C
#define HARRAY_C

#include "hArray.h"


///////////////////////////////////////////////////////
hArray::hArray(const char *name, TH1D* hseed, int nh){
  nameTag = name;
  hSeed = hseed;
  nHistos = nh;
  init();
}


///////////////////////////////////////////////////////
void hArray::init(){
  
  // initialize all from clones
  for (int i=0; i<nHistos; i++){
    TString hname = nameTag.Data();
    hname.Append(Form("_h%d",i));
    histos[i] = (TH1D*)hSeed->Clone(hname.Data());
    histos[i]->Reset();
  } 

  return; 
}



///////////////////////////////////////////////////////
void hArray::drawArray(){

  for (int ih=0; ih<nHistos; ih++){
    histos[ih]->SetLineColor(kCyan + 1);
  }

//  histos[0]->Draw("h");
  histos[0]->SetMinimum(0.);
  histos[0]->SetMaximum(1.3*histos[0]->GetMaximum());
  double binw = histos[0]->GetBinWidth(1);
  double buffer = 2.*binw;
  double xmin = histos[0]->GetBinLowEdge(1);
  double xmax = histos[0]->GetBinLowEdge(histos[0]->GetNbinsX())+binw;
  histos[0]->GetXaxis()->SetRangeUser(xmin+buffer,xmax-buffer);
  histos[0]->Draw("histo");
  for (int ih=1; ih<nHistos; ih++){
    histos[ih]->Draw("samehisto");
  }
  return;
}




///////////////////////////////////////////////////////
void hArray::delArray(){
  for (int ih=0; ih<nHistos; ih++){
    histos[ih]->Delete();
  }
  return;
}



///////////////////////////////////////////////////////
void hArray::setHistogram(int ih, TH1D* hh){
  TString hname = nameTag.Data();
  hname.Append(Form("_h%d",ih));
  histos[ih] = (TH1D*)hh->Clone(hname.Data());

  return;
}

#endif
