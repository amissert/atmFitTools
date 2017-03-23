#include "TH1D.h"
#include "TMath.h"
#include "TString.h"


using namespace std;


TH1D* GetCumulative(TH1D* h){

  // get name of new histogram
  TString hname = h->GetName();
  hname.Append("_CDF");

  // clone from original
  TH1D* hcum = (TH1D*)h->Clone(hname.Data());
  hcum->Reset();

  // fill with cum
  double sum = 0.;
  int N = h->GetNbinsX();
  for (int ibin=0; ibin<=h->GetNbinsX()+1; ibin++){
    hcum->SetBinContent(ibin,sum);
    sum += h->GetBinContent(ibin);
  }

  //
  return hcum;
}


TH1D* GetReverseCumulative(TH1D* h){

  // get name of new histogram
  TString hname = h->GetName();
  hname.Append("_CDF");

  // clone from original
  TH1D* hcum = (TH1D*)h->Clone(hname.Data());
  hcum->Reset();

  // fill with cum
  double sum = 0.;
  int N = h->GetNbinsX();
  for (int ibin=N+1; ibin>=0; ibin--){
    hcum->SetBinContent(ibin,sum);
    sum += h->GetBinContent(ibin);
  }

  //
  return hcum;
}


