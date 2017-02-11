#ifndef HBOUND_C
#define HBOUND_C

#include "TH1D.h"

/////////////////////////////////////////////////////////////
// Modify histogram filled by graph with physical lower bound
// Can be used as a template for enforcing an upper bound as well
void applyLoBound(TH1D* h, double lobound){
  
  // get bin width (assume constant);
  double binw = h->GetBinWidth(1);
  // identify the critical bin where the lower bound lies
  int critbin = -1;
  for (int ibin=1; ibin<=h->GetNbinsX(); ibin++){

     double binloval = h->GetBinLowEdge(ibin);
     double binhival  = binloval + h->GetBinWidth(ibin); 
     if ((binloval <= lobound) && (binhival> lobound)){
       critbin = ibin;
       break;
     }

  }

  // integrate un-physical bins
  double binsum=0;
  for (int ibin=(critbin-1); ibin>=0; ibin--){
    binsum+=h->GetBinContent(ibin);
  }
 
  // add to critical bin
  h->SetBinContent(critbin,h->GetBinContent(critbin)+binsum);

  // clear un-physical bins
  for (int ibin=(critbin-1); ibin>=0; ibin--){
    h->SetBinContent(ibin,0.);
  }
  
  //
  return;
}




#endif
