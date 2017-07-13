#ifndef  RCMap_C
#define  RCMap_C

#include "RCMap.h"


///////////////////////////////////////////////////
double RCMap::getRowSum(TH2D* h, int irow){
   int nbinsx = h->GetNbinsX();
   double sum = 0.;
   for (int xbin=1; xbin<=nbinsx; xbin++){
     int gbin = h->GetBin(xbin,irow);
     sum += h->GetBinContent(gbin);
   }
  
   return sum;
}


///////////////////////////////////////////////////
void RCMap::scaleRow(TH2D* h, int irow, double scale){
   int nbinsx = h->GetNbinsX();
   for (int xbin=1; xbin<=nbinsx; xbin++){
     int gbin = h->GetBin(xbin,irow);
     h->SetBinContent(gbin, (h->GetBinContent(gbin)*scale));
   }
   return ;
}



///////////////////////////////////////////////////
void RCMap::getConditionalX(){

   int nbinsx = hweight->GetNbinsX();
   int nbinsy = hweight->GetNbinsY();

   for (int ybin=1; ybin<=nbinsy; ybin++){
     double rowsum = getRowSum(hweight, ybin);
     double norm = 1.0;
     if (rowsum!=0){
       norm = 1./rowsum;
     }
     scaleRow(hweight,ybin,norm);
   }

   return;
}


///////////////////////////////////////////////////
double RCMap::get1RWeight(double rcpar){
  int ybin = hweight->GetYaxis()->FindBin(rcpar);
  return hweight->GetBinContent(2,ybin);
}


///////////////////////////////////////////////////
RCMap::RCMap(const char* fname){

  TFile* file = new TFile(fname);
  hmap = (TH2D*)file->Get("hmap");
  init();

}


void RCMap::init(){
  hweight = (TH2D*)hmap->Clone("hweight");
  getConditionalX();
}

#endif   /* ----- #ifndef RCMap_INC  ----- */

