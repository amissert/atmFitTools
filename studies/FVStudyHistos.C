
#include "TH1D.h"
#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TString.h"
#include "TGraph.h"
#include "TMath.h"
#include "THStack.h"
#include <iostream>
#include "../fqProcessedEvent.cxx"
#include "../eventSelectors.h"
#include "../TH2FV.h"

#ifndef NCAT
#define NCAT 6
#endif

#ifndef NSAMP 
#define NSAMP 2
#endif


//////////////////////////////////////////////////////////////////////////////////
// CLASS FVStudyHistos ///////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//  Container class for FV related histograms.  Filling histograms from processed
//  T2K events will probably become commonplace in the coming months.  This class
//  will organize these histograms and make for easy management.
//
//  Notes:
//    - Every histogram of interest has to indicies: [icat][isamp]
//        - "icat" is the event catagory. MC is broken down in to basic
//           catagories i.e. CCQE, CCnQE, CCMisID, ... 
//        - "isamp" is the event sample i.e. numu or nue
//    - "nameTag" is used to make unique histogram names
///////////////////////////////////////////////////////////////////////////////////
class FVStudyHistos{

  public:
  FVStudyHistos(const char* name);
  TString nameTag;


  TH1D* hwall[NCAT][NSAMP];
  TH1D* htowall[NCAT][NSAMP];
  TH1D* hwallv[NCAT][NSAMP];
  TH1D* htowallv[NCAT][NSAMP];
  TH1D* herec[NCAT][NSAMP];
  TH1D* hetrue[NCAT][NSAMP];
  TH1D* hmode[NCAT][NSAMP];
  TH1D* h[NCAT][NSAMP];

  TH2FV* hfv[NCAT][NSAMP];

  THStack* stack[2];

  double eventCount[NCAT][2];

  void init();
  void clearHistos();
  void drawStack(const char* var, int isample);
  void countEvent(int icat, int isamp, double ww);
  void printSummary(int isample);

};



///////////////////////////////////////////////////////////////////////////////////
void FVStudyHistos::printSummary(int isample){

  double sum=0; 
  sum += eventCount[0][isample];
  sum += eventCount[1][isample];
  sum += eventCount[2][isample];
  sum += eventCount[3][isample];
  sum += eventCount[4][isample];

  cout<<"-----------------------------"<<endl;
  cout<<"--EVENT SUMMARY--------------"<<endl;
  cout<<"-----------------------------"<<endl;
  cout<<"- CCQE:      "<<eventCount[0][isample]<<endl;
  cout<<"- CCnQE:     "<<eventCount[1][isample]<<endl;
  cout<<"- CCMisID:   "<<eventCount[2][isample]<<endl;
  cout<<"- NC:        "<<eventCount[3][isample]<<endl;
  cout<<"- Entering:  "<<eventCount[4][isample]<<endl;
  cout<<"----------------------------"<<endl;
  cout<<"- TOTAL:     "<<sum<<endl;
  cout<<endl;

  return;
}


///////////////////////////////////////////////////////////////////////////////////
void FVStudyHistos::countEvent(int icat, int isamp, double ww){

  eventCount[icat][isamp] += ww;

  return;
}



///////////////////////////////////////////////////////////////////////////////////
// draw all histograms in stack form with specified colors
void FVStudyHistos::drawStack(const char* var, int isample){

  // setup stack
  TString stackname = nameTag.Data();
  stackname.Append(Form("_stack%d",isample)); 
//  if (stack[isample]){
//    cout<<"FVStudyHistos: Deleting stack..."<<endl;
//    stack[isample]->Delete();
//    stack[isample] = new THStack(stackname.Data(),"");
//  }
//  else{
    cout<<"FVStudyHistos: Making new stack..."<<endl;
    stack[isample] = new THStack(stackname.Data(),"");
//  }

  // color scheme
  int colerz[NCAT] = {4,7,2,1,6};

  // variable name
  TString vname = var;

  // loop over catagories
  for (int icat=NCAT-1; icat>=0; icat--){

    TH1D* hadd;

    // get the right variabl
    if (!vname.CompareTo("wall")){
      hadd = hwall[icat][isample];
    }
    if (!vname.CompareTo("erec")){
      hadd = herec[icat][isample];
    }
    if (!vname.CompareTo("etrue")){
      hadd = hetrue[icat][isample];
    }
    if (!vname.CompareTo("wallv")){
      hadd = hwallv[icat][isample];
    }
    if (!vname.CompareTo("towall")){
      hadd = htowall[icat][isample];
    }
    if (!vname.CompareTo("towallv")){
      hadd = htowallv[icat][isample];
    }

    // make it the right color
    hadd->SetLineColor(colerz[icat]);
    hadd->SetFillColor(colerz[icat]);
    stack[isample]->Add(hadd);

  }

  // set some opts 
  stack[isample]->SetTitle(0);
  stack[isample]->Draw("h");

  //
  return;
}



///////////////////////////////////////////////////////////////////////////////////
// clear all histogram contents
void FVStudyHistos::clearHistos(){
  for (int isample=0; isample<NSAMP; isample++){
    for (int icat = 0; icat<NCAT; icat++){
       herec[icat][isample]->Reset();
       hetrue[icat][isample]->Reset();
       hwall[icat][isample]->Reset();
       hwallv[icat][isample]->Reset();
       htowall[icat][isample]->Reset();
       htowallv[icat][isample]->Reset();
    }
  }  
}


///////////////////////////////////////////////////////////////////////////////////
// initialize
void FVStudyHistos::init(){

   // name
  TString hname;
  for (int icat=0; icat<NCAT; icat++){
    for (int isamp=0; isamp<NSAMP; isamp++){
      hname = nameTag.Data();
      hname.Append(Form("h2fv_%d",isamp));
      hfv[icat][isamp] = new TH2FV(hname.Data(),-1,15,0,300,15,0,300);
    }
  }

  // for nbins
  int nbins = 40;

  // wall
  TH1D* wallseed = new TH1D("wallseed","wallseed",nbins,-100,500);
  wallseed->GetXaxis()->SetTitle("Wall [cm]");
  wallseed->SetStats(0);
  for (int icat=0; icat<NCAT; icat++){
    for (int isamp=0; isamp<NSAMP; isamp++){
      hname = nameTag.Data();
      hname.Append(Form("_hwall_samp%d_cat%d",isamp,icat));
      hwall[icat][isamp] = (TH1D*)wallseed->Clone(hname.Data());
      hwall[icat][isamp]->GetXaxis()->SetTitle("Wall [cm]");
      hwall[icat][isamp]->SetTitle(0);
      hwall[icat][isamp]->SetStats(0);
    }
  }

  // towall
  TH1D* towallseed = new TH1D("towallseed","towallseed",nbins,0,500);
  towallseed->GetXaxis()->SetTitle("Towall [cm]");
  for (int icat=0; icat<NCAT; icat++){
    for (int isamp=0; isamp<NSAMP; isamp++){
      hname = nameTag.Data();
      hname.Append(Form("_htowall_samp%d_cat%d",isamp,icat));
      htowall[icat][isamp] = (TH1D*)towallseed->Clone(hname.Data());
      htowall[icat][isamp]->GetXaxis()->SetTitle("Towall [cm]");
      htowall[icat][isamp]->SetTitle(0);
      htowall[icat][isamp]->SetStats(0);
    }
  }

  // wallv
  for (int icat=0; icat<NCAT; icat++){
    for (int isamp=0; isamp<NSAMP; isamp++){
      hname = nameTag.Data();
      hname.Append(Form("_hwallv_samp%d_cat%d",isamp,icat));
      hwallv[icat][isamp] = (TH1D*)wallseed->Clone(hname.Data());
      hwallv[icat][isamp]->GetXaxis()->SetTitle("True Wall [cm]");
      hwallv[icat][isamp]->SetTitle(0);
      hwallv[icat][isamp]->SetStats(0);
    }
  }

  // towallv
  for (int icat=0; icat<NCAT; icat++){
    for (int isamp=0; isamp<NSAMP; isamp++){
      hname = nameTag.Data();
      hname.Append(Form("_htowallv_samp%d_cat%d",isamp,icat));
      htowallv[icat][isamp] = (TH1D*)towallseed->Clone(hname.Data());
      htowallv[icat][isamp]->GetXaxis()->SetTitle("True Towall [cm]");
      htowallv[icat][isamp]->SetTitle(0);
      htowallv[icat][isamp]->SetStats(0);
    }
  }

  // erec
  TH1D* erecseed = new TH1D("erec","erec",nbins,0,2000);
  erecseed->GetXaxis()->SetTitle("E_{#nu}");
  for (int icat=0; icat<NCAT; icat++){
    for (int isamp=0; isamp<NSAMP; isamp++){
      hname = nameTag.Data();
      hname.Append(Form("=herec_samp%d_cat%d",isamp,icat));
      herec[icat][isamp] = (TH1D*)erecseed->Clone(hname.Data());
      herec[icat][isamp]->GetXaxis()->SetTitle("E_{#nu} [MeV]");
      herec[icat][isamp]->SetTitle(0);
      herec[icat][isamp]->SetStats(0);
    }
  }

  // etrue
  TH1D* etrueseed = new TH1D("etrue","etrue",nbins,0,2000);
  etrueseed->GetXaxis()->SetTitle("E_{#nu}");
  for (int icat=0; icat<NCAT; icat++){
    for (int isamp=0; isamp<NSAMP; isamp++){
      hname = nameTag.Data();
      hname.Append(Form("=hetrue_samp%d_cat%d",isamp,icat));
      hetrue[icat][isamp] = (TH1D*)etrueseed->Clone(hname.Data());
      hetrue[icat][isamp]->GetXaxis()->SetTitle("E_{#nu} [MeV]");
      hetrue[icat][isamp]->SetTitle(0);
      hetrue[icat][isamp]->SetStats(0);
    }
  }

  // event count
  for (int icat=0; icat<NCAT; icat++){
    for (int isamp=0; isamp<NSAMP; isamp++){
      eventCount[icat][isamp] = 0.;
    }
  }

  //
  return;
}



///////////////////////////////////////////////////////////////////////////////////
// construct
FVStudyHistos::FVStudyHistos(const char* name){
 nameTag = name;
 init();
}




















