
#include "TH1D.h"


TH1D* getfracerr(TH1D* h){

  TH1D* herr = (TH1D*)h->Clone("herr");

  herr->Reset();

  for (int i=0; i<h->GetNbinsX(); i++){
    double binc = h->GetBinContent(i);
    double binerr = h->GetBinError(i);
    if (binerr!=0){
      herr->SetBinContent(i,(binerr*binerr)/binc);
    }
  }

  return herr;
}
