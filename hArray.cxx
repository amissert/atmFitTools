
#ifndef HARRAY_C
#define HARRAY_C

#include "hArray.h"


///////////////////////////////////////////////////////
void hArray::calcSummary(){

  // get template
  TString hname = nameTag.Data();
  hname.Append("_array_summary");
  hSummary = (TH1D*)hSeed->Clone(hname.Data());
  hSummary->Reset();

  // calc mean and rms for each bin
  const int N = nHistos;
  double contents[N];
  for (int ibin=1; ibin<=hSummary->GetNbinsX(); ibin++){
    
    // fill array
    for (int ipt=0; ipt<N; ipt++){
      contents[ipt] = histos[ipt]->GetBinContent(ibin);
    }

    // get array stats
    double mean = arraymeanD(contents,N);
    double sig  = TMath::Sqrt(arrayvarD(contents,N))/1.;

    // set as summary
    hSummary->SetBinContent(ibin,mean);
    hSummary->SetBinError(ibin,sig);
  }

  hSummary->SetFillColor(kCyan + 1);
  hSummary->SetMarkerStyle(1);
  double binw = hSummary->GetBinWidth(1);
  double buffer = 3.*binw;
  double xmin = hSummary->GetBinLowEdge(1);
  double xmax = hSummary->GetBinLowEdge(hSummary->GetNbinsX())+binw;
  hSummary->GetXaxis()->SetRangeUser(xmin+buffer,xmax-buffer);
  return;
}

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
  double buffer = 4.*binw;
  double xmin = histos[0]->GetBinLowEdge(1);
  double xmax = histos[0]->GetBinLowEdge(histos[0]->GetNbinsX())+binw;
  histos[0]->GetXaxis()->SetRangeUser(xmin+buffer,xmax-buffer);
  histos[0]->Draw("histo");
  for (int ih=1; ih<nHistos; ih++){
    histos[ih]->SetLineWidth(1);
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
