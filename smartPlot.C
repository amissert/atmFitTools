#ifndef SMARTPLOT_C
#define SMARTPLOT_C

#include "TH2D.h"
#include "TH2Poly.h"
#include "TH2FV.h"


///////////////////////////////////////////////////////////////////
void smartAxis(TH2FV* h){
 
  double min=1e11;
  double max=-1e11;
  double buffsize = 0.1; 
  for (int ibin =0; ibin<h->GetNumberOfBins(); ibin++){
     double content = h->GetBinContent(ibin);
     if (content ==0) continue;
     if (content < min){
       min = content;
     }
     if (content > max){
       max = content;
     }
  }

  double range = max - min;
  double buffer = range*buffsize;
  h->SetMinimum(min - buffer);
  h->SetMaximum(max + buffer);

  //
  return;
}



#endif

