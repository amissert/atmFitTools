#include "TH2D.h"
#include "TString.h"
#include "TRandom2.h"
#include <iostream>

using namespace std;

class softcompare{
  public:

  softcompare();

  TH2D* makemask(TH2D*  hh, double mincontent);
  TH2D* softratio(TH2D* hmc, TH2D* hdata, double zmin=0, TH2D* hstats=NULL);
  TH2D* softdiff(TH2D* hmc, TH2D* hdata, double zmin=0);
//  TH2D* softfrac(TH2D* hmc, TH2D* hdata, double zmin=15);

  TH2D* htest[2];
  TH2D* htestratio;
  void test(int nev,double minc);

};

void softcompare::test(int nev, double minc){

  TRandom2* randy = new TRandom2(nev);

  htest[0] = new TH2D("h1","h1",50,0,10,50,0,10);
  htest[1] = new TH2D("h2","h2",50,0,10,50,0,10);

  for (int i=0; i<nev;  i++){
    htest[0]->Fill(randy->Gaus(4,1),randy->Gaus(4,1));
    htest[1]->Fill(randy->Gaus(6,1),randy->Gaus(6,1));
  }

  htestratio = softdiff(htest[0],htest[1],minc);

  htestratio->Draw("colz");

}

// mask out all bins with content below some threshold
TH2D* softcompare::makemask(TH2D*  hh, double mincontent){

  // clone a new histogram
  TString hname = hh->GetName();
  hname.Append("_mask");
  TH2D* hmask = (TH2D*)hh->Clone(hname.Data());
  hmask->Reset();

  // fill mask
  for (int xbin=1; xbin<=hmask->GetNbinsX(); xbin++){
    for (int ybin=1; ybin<=hmask->GetNbinsY(); ybin++){
       double content = hh->GetBinContent(xbin,ybin);
       if (content<mincontent) hmask->SetBinContent(xbin,ybin,0.);
       else{
         hmask->SetBinContent(xbin,ybin,1.);
       }
    }
  }

  // return mask
  return hmask;
}


TH2D* softcompare::softdiff(TH2D* hmc, TH2D* hdata, double zmin){

  // clone a new histogram
  TString hname = hmc->GetName();
  hname.Append("_softratio");
  TH2D* hratio = (TH2D*)hmc->Clone(hname.Data());
  hratio->Reset();

  // make masks
  TH2D* hmasktot = makemask(hmc,zmin);
  TH2D* hmaskdata = makemask(hdata,zmin);
  hmasktot->Multiply(hmaskdata);

  // make it soft
  hmc->Smooth();
  hdata->Smooth();

  // scale it
  double norm = hmc->Integral()/hdata->Integral();
  hmc->Scale(1./norm);

  // subtract it
  hdata->Scale(-1.);
  hmc->Add(hdata);

  // mask it
  hmc->Multiply(hmasktot);
  hmc->SetContour(100);

  // return it
  return hmc;

}

TH2D* softcompare::softratio(TH2D* h1, TH2D* h2, double zmin, TH2D* hstats){

  // clone a new histogram
  TString hname = h1->GetName();
  hname.Append("_softratio");
  TH2D* hratio = (TH2D*)h1->Clone(hname.Data());
  hratio->Reset();

  // clone inputs
  TH2D* hnum = (TH2D*)h1->Clone("hnum"); 
  TH2D* hdenom = (TH2D*)h2->Clone("hdenom"); 

  // make it soft
//  hnum->Smooth();
  hnum->Smooth();
//  hdenom->Smooth();
  hdenom->Smooth();

  // scale it
  double norm = hnum->Integral()/hdenom->Integral();
  cout<<"norm: "<<norm<<endl;
  hnum->Scale(1./norm);
  hdenom->Add(hnum);

  // mask low statistics
  TH2D* hmask;
  if (hstats==NULL) hmask = makemask(hdenom,zmin);
  else{
    hmask = makemask(hstats,zmin);
  }

  hdenom->Multiply(hmask);
  hnum->Multiply(hmask);


  // divide it
  hnum->Divide(hdenom);

  // mask it
  hnum->SetContour(100);

  // return it
  hnum->Draw("colz");
  return hnum;

}

softcompare::softcompare(){


}
