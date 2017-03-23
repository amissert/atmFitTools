#ifndef FVSTUD_C
#define FVSTUD_C

#include "TH1D.h"
#include "TH2D.h"
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
#define NSAMP 3
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

  // histograms
  TH1D* hwall[NCAT][NSAMP];
  TH1D* htowall[NCAT][NSAMP];
  TH1D* hwallv[NCAT][NSAMP];
  TH1D* htowallv[NCAT][NSAMP];
  TH1D* herec[NCAT][NSAMP];
  TH1D* hetrue[NCAT][NSAMP];
  TH1D* hmode[NCAT][NSAMP];
  TH1D* h[NCAT][NSAMP];
  TH1D* heres[NSAMP];
  TH1D* hdirres[NSAMP];
  TH2D* heresvwall[NSAMP]; 
  TH2D* heresvtowall[NSAMP]; 
  TH2D* hpos[NCAT][NSAMP]; 
  TH1D* herectot[NSAMP];
  TH2FV* hfv[NCAT][NSAMP];
  TH2FV* hfvtot[NSAMP];
  TH2FV* hfvtottrue[NSAMP];
  TH2FV* hfveres[NSAMP];
  TH2FV* hfveresfrac[NSAMP];
  THStack* stack[3];
  double eventCount[NCAT][3];

  void init();
  void clearHistos();
  void drawStack(const char* var, int isample);
  void countEvent(int icat, int isamp, double ww);
  void printSummary(int isample);
  void normalizeResiduals();

};


///////////////////////////////////////////////////////////////////////////////////
void FVStudyHistos::normalizeResiduals(){

  double norm = 1.0;
  norm = 1./heres[0]->Integral();
  heres[0]->Scale(norm);
  //
  norm = 1./heres[1]->Integral();
  heres[1]->Scale(norm);
  //
  norm = 1./hdirres[0]->Integral();
  hdirres[0]->Scale(norm);
  //
  norm = 1./hdirres[1]->Integral();
  hdirres[2]->Scale(norm);

  return;
}

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
  cout<<"FVStudyHistos: Making new stack..."<<endl;
  stack[isample] = new THStack(stackname.Data(),"");

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
  stack[isample]->Draw("histo");

  //
  if (!vname.CompareTo("towallv")){
    stack[isample]->GetXaxis()->SetTitle("True Towall [cm]");
  }
  if (!vname.CompareTo("wallv")){
    stack[isample]->GetXaxis()->SetTitle("True Wall [cm]");
  }
  if (!vname.CompareTo("towall")){
    stack[isample]->GetXaxis()->SetTitle("Towall [cm]");
  }
  if (!vname.CompareTo("wall")){
    stack[isample]->GetXaxis()->SetTitle("Wall [cm]");
  }
  if (!vname.CompareTo("erec")){
    stack[isample]->GetXaxis()->SetTitle("E_{#nu} [MeV]");
  }
  if (!vname.CompareTo("etrue")){
    stack[isample]->GetXaxis()->SetTitle("True E_{#nu} [MeV]");
  }

  stack[isample]->Draw("histo");

  //
  return;
}



///////////////////////////////////////////////////////////////////////////////////
// clear all histogram contents
void FVStudyHistos::clearHistos(){

  // clear the histos
  for (int isample=0; isample<NSAMP; isample++){
    for (int icat = 0; icat<NCAT; icat++){
       herec[icat][isample]->Reset();
       hetrue[icat][isample]->Reset();
       hwall[icat][isample]->Reset();
       hwallv[icat][isample]->Reset();
       htowall[icat][isample]->Reset();
       htowallv[icat][isample]->Reset();
       hpos[icat][isample]->Reset();
       // clear the count too
       eventCount[icat][isample] = 0.;
    }
    heres[isample]->Reset();
    hdirres[isample]->Reset();
  }  

  return;
}


///////////////////////////////////////////////////////////////////////////////////
// initialize
void FVStudyHistos::init(){

  // name
  TString hname;
  int nfvbins = 30;

  // 2D histograms 
  for (int icat=0; icat<NCAT; icat++){
    for (int isamp=0; isamp<NSAMP; isamp++){
      hname = nameTag.Data();
      hname.Append(Form("h2fv_sample%d_cat%d",isamp,icat));
      hfv[icat][isamp] = new TH2FV(hname.Data(),-1,nfvbins,0,300,nfvbins,0,300);
      hfv[icat][isamp]->SetName(hname.Data());
      // position
      hname = nameTag.Data();
      hname.Append(Form("hposition_samp%d_cat%d",isamp,icat));
      hpos[icat][isamp] = new TH2D(hname.Data(),hname.Data(),200,0,1800,200,0,2200);
      hpos[icat][isamp]->SetName(hname.Data());
    }
  }
  
  // hfv totals
  hname = nameTag.Data();
  hname.Append("hfvtot_nue");
  hfvtot[0] = new TH2FV(hname.Data(),-1,nfvbins,0,300,nfvbins,0,300);
  hfvtot[0]->SetName(hname.Data());
  hname = nameTag.Data();
  hname.Append("hfvtot_numu");
  hfvtot[1] = new TH2FV(hname.Data(),-1,nfvbins,0,300,nfvbins,0,300);
  hfvtot[1]->SetName(hname.Data());
  hname = nameTag.Data();
  hname.Append("hfvtot_nue1rpi");
  hfvtot[2] = new TH2FV(hname.Data(),-1,nfvbins,0,300,nfvbins,0,300);
  hfvtot[2]->SetName(hname.Data());
  //
  hname = nameTag.Data();
  hname.Append("hfvtottrue_nue");
  hfvtottrue[0] = new TH2FV(hname.Data(),-1,nfvbins,0,300,nfvbins,0,300);
  hfvtottrue[0]->SetName(hname.Data());
  hname = nameTag.Data();
  hname.Append("hfvtottrue_numu");
  hfvtottrue[1] = new TH2FV(hname.Data(),-1,nfvbins,0,300,nfvbins,0,300);
  hfvtottrue[1]->SetName(hname.Data());
  hname = nameTag.Data();
  hname.Append("hfvtottrue_nue1rpi");
  hfvtottrue[2] = new TH2FV(hname.Data(),-1,nfvbins,0,300,nfvbins,0,300);
  hfvtottrue[2]->SetName(hname.Data());

  // hfv eres 
  hname = nameTag.Data();
  hname.Append("hfveres_nue");
  hfveres[0] = new TH2FV(hname.Data(),-1,nfvbins,0,300,nfvbins,0,300);
  hfveres[0]->SetName(hname.Data());
  hname = nameTag.Data();
  hname.Append("hfveres_numu");
  hfveres[1] = new TH2FV(hname.Data(),-1,nfvbins,0,300,nfvbins,0,300);
  hfveres[1]->SetName(hname.Data());
  hname = nameTag.Data();
  hname.Append("hfveres_nue1rpi");
  hfveres[2] = new TH2FV(hname.Data(),-1,nfvbins,0,300,nfvbins,0,300);
  hfveres[2]->SetName(hname.Data());

  // hfv eres fractional
  hname = nameTag.Data();
  hname.Append("hfveresfrac_nue");
  hfveresfrac[0] = new TH2FV(hname.Data(),-1,nfvbins,0,300,nfvbins,0,300);
  hfveresfrac[0]->SetName(hname.Data());
  hname = nameTag.Data();
  hname.Append("hfveresfrac_numu");
  hfveresfrac[1] = new TH2FV(hname.Data(),-1,nfvbins,0,300,nfvbins,0,300);
  hfveresfrac[1]->SetName(hname.Data());
  hname = nameTag.Data();
  hname.Append("hfveresfrac_nue1rpi");
  hfveresfrac[2] = new TH2FV(hname.Data(),-1,nfvbins,0,300,nfvbins,0,300);
  hfveresfrac[2]->SetName(hname.Data());

  //
  hname = nameTag.Data();
  hname.Append("hfvtottrue_nue");
  hfvtottrue[0] = new TH2FV(hname.Data(),-1,nfvbins,0,300,nfvbins,0,300);
  hfvtottrue[0]->SetName(hname.Data());
  hname = nameTag.Data();
  hname.Append("hfvtottrue_numu");
  hfvtottrue[1] = new TH2FV(hname.Data(),-1,nfvbins,0,300,nfvbins,0,300);
  hfvtottrue[1]->SetName(hname.Data());
  hname = nameTag.Data();
  hname.Append("hfvtottrue_nue1rpi");
  hfvtottrue[2] = new TH2FV(hname.Data(),-1,nfvbins,0,300,nfvbins,0,300);
  hfvtottrue[2]->SetName(hname.Data());

  // for nbins
  int nbins = 40;
  double wallmax = 500.;
  // wall
  TH1D* wallseed = new TH1D("wallseed","wallseed",nbins,-100,wallmax);
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
  TH1D* towallseed = new TH1D("towallseed","towallseed",nbins,0,wallmax);
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
  double maxE[2] = {1500.,2500.};
  TH1D* erecseed[3];
  erecseed[0] = new TH1D("erecnue","erecnue",nbins,0,maxE[0]);
  erecseed[1] = new TH1D("erecnumu","erecnumu",nbins,0,maxE[1]);
  erecseed[2] = new TH1D("erecnue1rpi","erecnue1rpi",nbins/2,0,maxE[0]);
  erecseed[0]->GetXaxis()->SetTitle("E_{#nu}");
  erecseed[1]->GetXaxis()->SetTitle("E_{#nu}");
  erecseed[2]->GetXaxis()->SetTitle("E_{#nu}");
  for (int icat=0; icat<NCAT; icat++){
    for (int isamp=0; isamp<NSAMP; isamp++){
      hname = nameTag.Data();
      hname.Append(Form("_herec_samp%d_cat%d",isamp,icat));
      herec[icat][isamp] = (TH1D*)erecseed[isamp]->Clone(hname.Data());
      herec[icat][isamp]->GetXaxis()->SetTitle("E_{#nu} [MeV]");
      herec[icat][isamp]->SetTitle(0);
      herec[icat][isamp]->SetStats(0);
    }
  }
  herectot[0] = new TH1D("erece","erece",nbins,0,maxE[0]);
  herectot[1] = new TH1D("erecmu","erecmu",nbins,0,maxE[1]);
  herectot[2] = new TH1D("erec1rpi","erece1rpi",nbins,0,maxE[1]);


  // etrue
  for (int icat=0; icat<NCAT; icat++){
    for (int isamp=0; isamp<NSAMP; isamp++){
      hname = nameTag.Data();
      hname.Append(Form("_hetrue_samp%d_cat%d",isamp,icat));
      hetrue[icat][isamp] = (TH1D*)erecseed[isamp]->Clone(hname.Data());
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

  // residuals
  for (int isamp=0; isamp<NSAMP; isamp++){
    hname = nameTag.Data();
    hname.Append(Form("_heresidual_samp%d",isamp));
    heres[isamp]  = new TH1D(hname.Data(),hname.Data(),80,-1000,1000);
    hname = nameTag.Data();
    hname.Append(Form("_hdirresidual_samp%d",isamp));
    hdirres[isamp]  = new TH1D(hname.Data(),hname.Data(),200,0,0.5);
    hname = nameTag.Data();
    hname.Append(Form("_heres_v_wall_samp%d",isamp));
    heresvwall[isamp] = new TH2D(hname.Data(),hname.Data(),60,0,300,100,-1000,1000);
    hname = nameTag.Data();
    hname.Append(Form("_heres_v_towall_samp%d",isamp));
    heresvtowall[isamp] = new TH2D(hname.Data(),hname.Data(),60,0,300,100,-1000,1000);
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



















#endif
