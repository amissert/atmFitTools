#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TChain.h";
#include "TTree.h"
#include "TLine.h"
#include <iostream>
#include "TLegend.h"
#include "TRandom2.h"


using namespace std;

class markovCompare{
  public:

  markovCompare(const char* fname1 , const char* fname2);

  TTree* partree1;
  TTree* partree2;
  TFile* parfile1;
  TFile* parfile2;
  TH1D* h1;
  TH1D* h2;
  TH2D* h2a;
  TH2D* h2b;
  TLine* tline;
  TLine* tline2;
  TCanvas* cc;
  int nburn1;
  int nburn2;

  void compare1DPar(int ipar,double min, double max, double partrue);
  void compare2DPar(int ipar1, double xmin1, double xmax1, double partrue1,
                                 int ipar2, double xmin2, double xmax2, double partrue2);
};


void markovCompare::compare2DPar(int ipar1, double xmin1, double xmax1, double partrue1,
                                 int ipar2, double xmin2, double xmax2, double partrue2){

 // make histograms
 h2a = new TH2D("h2a","h2a",100,xmin1,xmax1,100,xmin2,xmax2); 
 h2b = new TH2D("h2b","h2b",100,xmin1,xmax1,100,xmin2,xmax2); 

 // set par array and set branch address
 double par[5000];
 partree1->SetBranchAddress("par",par);
 for (int i=nburn1; i<partree1->GetEntries(); i++){
   partree1->GetEntry(i);
   h2a->Fill(par[ipar1],par[ipar2]);
 }
 partree2->SetBranchAddress("par",par);
 for (int i=nburn2; i<partree2->GetEntries(); i++){
   partree2->GetEntry(i);
   h2b->Fill(par[ipar1],par[ipar2]);
 }

 // True value lines
 tline = new TLine(partrue1,0,partrue1,100000); 
 tline->SetLineColor(kOrange);
 tline->SetLineWidth(3);
 tline2 = new TLine(0,partrue2,100000,partrue2); 
 tline2->SetLineColor(kOrange);
 tline2->SetLineWidth(3);
// h2->SetLineColor(kBlue);
 h2a->Scale(h2b->Integral()/h2a->Integral());

 h2b->SetFillColor(kBlue);
 h2b->SetLineColor(kBlue);
 h2a->SetLineColor(kBlack);
 h2a->SetFillColor(kBlack);
 h2b->Draw("box");
 h2a->Draw("boxsame");
// if (h1->GetMaximum()>h2->GetMaximum()){
//   h1->Draw();
//   h2->Draw("same");
// }
// else{
//   h2->Draw();
 //  h1->Draw("same");
// }
// tline->Draw("same");
 return;
}

void markovCompare::compare1DPar(int ipar, double xmin, double xmax, double partrue){

 int thinning = 10; 

 h1 = new TH1D("h1","h1",100,xmin,xmax); 
 h2 = new TH1D("h2","h2",100,xmin,xmax); 


 double par[5000];
 partree1->SetBranchAddress("par",par);
 for (int i=nburn1; i<partree1->GetEntries(); i++){
   if ((i%thinning)!=0) continue;
   partree1->GetEntry(i);
   h1->Fill(par[ipar]);
 }


 partree2->SetBranchAddress("par",par);
 for (int i=nburn2; i<partree2->GetEntries(); i++){
   if ((i%thinning)!=0) continue;
   partree2->GetEntry(i);
   h2->Fill(par[ipar]);
 }

 tline = new TLine(partrue,0,partrue,100000); 
 tline->SetLineColor(kOrange);
 tline->SetLineWidth(3);
 h2->SetLineColor(kBlue);
 h2->Scale(h1->Integral()/h2->Integral());

 if (h1->GetMaximum()>h2->GetMaximum()){
   h1->Draw();
   h2->Draw("same");
 }
 else{
   h2->Draw();
   h1->Draw("same");
 }
 tline->Draw("same");
 return;
}

markovCompare::markovCompare(const char* fname1, const char* fname2){

 // cc = new TCanvas("cc","cc",800,700);

  parfile1 = new TFile(fname1);
  parfile2 = new TFile(fname2);

  partree1 = (TTree*)parfile1->Get("MCMCpath");
  partree2 = (TTree*)parfile2->Get("MCMCpath");
  
}
