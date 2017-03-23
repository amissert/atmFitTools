#ifndef T2KFV_C
#define T2KFV_C

#include "TH1D.h"
#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TGraph.h"
#include "TMath.h"
#include "THStack.h"
#include "TVector3.h"
#include "TObjArray.h"
#include "TCanvas.h"
#include <iostream>
#include "../fqProcessedEvent.cxx"
#include "../eventSelectors.h"
#include "../TH2FV.h"
#include "FVStudyHistos.C"
#include "GetCumulative.C"

#ifndef NVOLS 
#define NVOLS 1
#endif





//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
class T2KFVStudies{

  public:

  T2KFVStudies(TChain* ch1);
  TChain *chMC;
  fqProcessedEvent* mcEvent;
  fqcutparams fqpars;

  FVStudyHistos* histos[NVOLS]; //< histos for each FV volume

  FVStudyHistos* histosCompare[5]; //< pass/fail histos for comparing cuts 

  TCanvas* canvas;

  THStack* stack[10];
  TH2D* hpoorfv[2];
  TH1D*  hpoor[2];
  TH1D*  hpoorlmom[2];
  TH1D*  hpoortw[2];

  int Nmax;
//  int Nmax;
  
  // code for fv cuts
  // 0 -> no wall or towal cuts
  // 1 -> nominal wall cuts
  // 2 -> inverse nominal cuts
  // 3 -> small wall cut
  // 4 -> entering bg cut
//  int fvcode;
  int passedCuts;
  int passedMuCuts;
  int passedECuts;
  int passedE1RPiCuts;
  double NormPOT;

  // loop
  void loop();
  int passT2KCuts();
  int passFVCuts(int ivol);

  // fill a particular instance of FVStudyHistos
  void fill(FVStudyHistos* hh);
  void fillHistos(int ivol);
  void compareCuts(double towall1, double wall1, double towall2, double wall2);
  int getCatagory(int wantnutype);
  void drawSummary(int ivol, int isample);
  void drawFVSummary(int icat, int ivol, int isample);
  void drawCompare(const char* var);
  void drawPassFail(const char* var, int isamp);
  void findPoorTowall();
  float calcResPiNuEnergy(float lmom, const TVector3* dir);

  // print many histograms to a directory
  void printComparePlots(const char* directory, const char* tag = "compare_");
  void printFVPlots(const char* directory);
  // histogram manipulation tools
  double stackIntegral(THStack* hstack);
  void makeAxesSame(THStack* stack1, THStack* stack2);
  void init();
};



//////////////////////////////////////////////////////////////////////////////////////////
float T2KFVStudies::calcResPiNuEnergy(float lmom, const TVector3* dir){

  // values
  float Mp = 938.272;
  float V  = 27.0;
  float MDelta = 1232.0;
  float Me = 0.511;
  float Ee = TMath::Sqrt(Me*Me + lmom*lmom);

  // cosine to beam
  float costh = 0.669764*dir->X() + -0.742179*dir->Y() + 0.024228*dir->Z();

  // RC energy
  float Erec = (2.*Mp*Ee + MDelta*MDelta - Mp*Mp - Me*Me)
              /(2.*(Mp - Ee + lmom*costh));

  //
  return Erec;
}



//////////////////////////////////////////////////////////////////////////////////////////
void T2KFVStudies::printFVPlots(const char* dir){

  TString plotdir = dir;
  TString pname = "";
  
  TCanvas* cc = new TCanvas("cc","cc",800,700);

  // total fv density
  pname = plotdir.Data();
  pname.Append("hfv_nue_total.png");
  histos[0]->hfvtot[0]->Draw("colz");
  cc->Print(pname.Data());
  //
  pname = plotdir.Data();
  pname.Append("hfv_numu_total.png");
  histos[0]->hfvtot[1]->Draw("colz");
  cc->Print(pname.Data());
  //

  // catagories
  for (int icat=0; icat<NCAT; icat++){
   
    // fv
    pname = plotdir.Data();
    pname.Append(Form("hfv_nue_cat%d.png",icat));
    histos[0]->hfv[icat][0]->Draw("colz");
    cc->Print(pname.Data());
    //
    pname = plotdir.Data();
    pname.Append(Form("hfv_numu_cat%d.png",icat));
    histos[0]->hfv[icat][1]->Draw("colz");
    cc->Print(pname.Data());
    //

    // pos
    pname = plotdir.Data();
    pname.Append(Form("hpos_nue_cat%d.png",icat));
    histos[0]->hpos[icat][0]->Draw("colz");
    cc->Print(pname.Data());
    //
    pname = plotdir.Data();
    pname.Append(Form("hpos_numu_cat%d.png",icat));
    histos[0]->hpos[icat][1]->Draw("colz");
    cc->Print(pname.Data());
    //
  }
  pname = plotdir.Data();
  pname.Append("hwall_numu.png");
  histos[0]->drawStack("wall",1);
  cc->Print(pname.Data());
  //
  pname = plotdir.Data();
  pname.Append("htowall_numu.png");
  histos[0]->drawStack("towall",1);
  cc->Print(pname.Data());
  //
  pname = plotdir.Data();
  pname.Append("hwall_nue.png");
  histos[0]->drawStack("wall",0);
  cc->Print(pname.Data());
  //
  pname = plotdir.Data();
  pname.Append("htowall_nue.png");
  histos[0]->drawStack("towall",0);
  cc->Print(pname.Data());
  
  //
  return;
}



//////////////////////////////////////////////////////////////////////////////////////////
void T2KFVStudies::printComparePlots(const char* directory, const char* tag){
  
  TString stag = tag;
  TString sdir = directory;

  TString pltname;

  // draw pass vs fail plots
  pltname = sdir.Data(); 
  pltname.Append(stag.Data());
  pltname.Append("passfail_erec_nue.png");
  drawPassFail("erec",0);
  canvas->Print(pltname.Data());
  //
  pltname = sdir.Data(); 
  pltname.Append(stag.Data());
  pltname.Append("passfail_erec_numu.png");
  drawPassFail("erec",1);
  canvas->Print(pltname.Data());
  //
  pltname = sdir.Data(); 
  pltname.Append(stag.Data());
  pltname.Append("passfail_erec_nue1rpi.png");
  drawPassFail("erec",2);
  canvas->Print(pltname.Data());


  pltname = sdir.Data(); 
  pltname.Append(stag.Data());
  pltname.Append("passfail_wall_nue.png");
  drawPassFail("wall",0);
  canvas->Print(pltname.Data());
  //
  pltname = sdir.Data(); 
  pltname.Append(stag.Data());
  pltname.Append("passfail_wall_numu.png");
  drawPassFail("wall",1);
  canvas->Print(pltname.Data());
  //
  pltname = sdir.Data(); 
  pltname.Append(stag.Data());
  pltname.Append("passfail_wall_nue1rpi.png");
  drawPassFail("wall",2);
  canvas->Print(pltname.Data());


  pltname = sdir.Data(); 
  pltname.Append(stag.Data());
  pltname.Append("passfail_wallv_nue.png");
  drawPassFail("wallv",0);
  canvas->Print(pltname.Data());
  //
  pltname = sdir.Data(); 
  pltname.Append(stag.Data());
  pltname.Append("passfail_wallv_numu.png");
  drawPassFail("wallv",1);
  canvas->Print(pltname.Data());
  //
  pltname = sdir.Data(); 
  pltname.Append(stag.Data());
  pltname.Append("passfail_wallv_nue1rpi.png");
  drawPassFail("wallv",2);
  canvas->Print(pltname.Data());


  //
  pltname = sdir.Data(); 
  pltname.Append(stag.Data());
  pltname.Append("passfail_etrue_nue.png");
  drawPassFail("etrue",0);
  canvas->Print(pltname.Data());
  //
  pltname = sdir.Data(); 
  pltname.Append(stag.Data());
  pltname.Append("passfail_etrue_numu.png");
  drawPassFail("etrue",1);
  canvas->Print(pltname.Data());
  //
  pltname = sdir.Data(); 
  pltname.Append(stag.Data());
  pltname.Append("passfail_etrue_nuepi.png");
  drawPassFail("etrue",2);
  canvas->Print(pltname.Data());


  pltname = sdir.Data(); 
  pltname.Append(stag.Data());
  pltname.Append("passfail_towall_nue.png");
  drawPassFail("towall",0);
  canvas->Print(pltname.Data());
  //
  pltname = sdir.Data(); 
  pltname.Append(stag.Data());
  pltname.Append("passfail_towall_numu.png");
  drawPassFail("towall",1);
  canvas->Print(pltname.Data());
  //
  pltname = sdir.Data(); 
  pltname.Append(stag.Data());
  pltname.Append("passfail_towall_nue1rpi.png");
  drawPassFail("towall",2);
  canvas->Print(pltname.Data());


  pltname = sdir.Data(); 
  pltname.Append(stag.Data());
  pltname.Append("passfail_towallv_nue.png");
  drawPassFail("towallv",0);
  canvas->Print(pltname.Data());
  //
  pltname = sdir.Data(); 
  pltname.Append(stag.Data());
  pltname.Append("passfail_towallv_numu.png");
  drawPassFail("towallv",1);
  canvas->Print(pltname.Data());
  //
  pltname = sdir.Data(); 
  pltname.Append(stag.Data());
  pltname.Append("passfail_towallv_nue.png");
  drawPassFail("towallv",2);
  canvas->Print(pltname.Data());


  // residual in difference
  canvas = new TCanvas("canvas","canvas",800,700);
  pltname = sdir.Data(); 
  pltname.Append(stag.Data());
  pltname.Append("difference_erec_nue.png");
  histosCompare[4]->drawStack("erec",0);
  canvas->Print(pltname.Data());
  //
  pltname = sdir.Data(); 
  pltname.Append(stag.Data());
  pltname.Append("difference_erec_numu.png");
  histosCompare[4]->drawStack("erec",1);
  canvas->Print(pltname.Data());
  //
  pltname = sdir.Data(); 
  pltname.Append(stag.Data());
  pltname.Append("difference_erec_nue.png");
  histosCompare[4]->drawStack("erec",2);
  canvas->Print(pltname.Data());

  //
  pltname = sdir.Data(); 
  pltname.Append(stag.Data());
  pltname.Append("difference_eres_nue.png");
  histosCompare[4]->heres[0]->SetLineColor(kRed);
  histosCompare[4]->heres[0]->SetLineWidth(3);
  histosCompare[4]->heres[0]->Draw("histo");
  histosCompare[0]->heres[0]->SetLineColor(kBlack);
  histosCompare[0]->heres[0]->SetLineWidth(3);
  histosCompare[0]->heres[0]->Draw("samehist");
  canvas->Print(pltname.Data());
  //
  pltname = sdir.Data(); 
  pltname.Append(stag.Data());
  pltname.Append("difference_eres_numu.png");
  histosCompare[4]->heres[1]->SetLineColor(kRed);
  histosCompare[4]->heres[1]->SetLineWidth(3);
  histosCompare[4]->heres[1]->Draw("h");
  histosCompare[0]->heres[1]->SetLineColor(kBlack);
  histosCompare[0]->heres[1]->SetLineWidth(3);
  histosCompare[0]->heres[1]->Draw("samehisto");
  canvas->Print(pltname.Data());
  //
  pltname = sdir.Data(); 
  pltname.Append(stag.Data());
  pltname.Append("difference_eres_nue1rpi.png");
  histosCompare[4]->heres[2]->SetLineColor(kRed);
  histosCompare[4]->heres[2]->SetLineWidth(3);
  histosCompare[4]->heres[2]->Draw("h");
  histosCompare[0]->heres[2]->SetLineColor(kBlack);
  histosCompare[0]->heres[2]->SetLineWidth(3);
  histosCompare[0]->heres[2]->Draw("samehisto");
  canvas->Print(pltname.Data());
  
  return;

}




//////////////////////////////////////////////////////////////////////////////////////////
void T2KFVStudies::makeAxesSame(THStack* stack1, THStack* stack2){

  // get the maximums
  double max1 = stack1->GetMaximum();
  double max2 = stack2->GetMaximum();

  // find the highests max
  double gmax = TMath::Max(max1,max2);

  // set ranges
  stack1->SetMinimum(0);
  stack1->SetMaximum(gmax);
  //
  stack2->SetMinimum(0);
  stack2->SetMaximum(gmax);
 
  return;
}



/*
//////////////////////////////////////////////////////////////////////////////////////////
void T2KFVStudies::drawStack(const char* var, int istack){

  // setup stack
  stack[istack] = new THStack(Form("hs%d",is),Form("hstack%d",is));

  // color scheme
  int colerz[NCAT] = {4,7,2,1,6};

  // loop over catagories
  for (int icat=NCAT-1; icat>=0; icat--){
     histos[ivol]->hwall[icat][isample]->SetLineColor(colerz[icat]);
     histos[ivol]->hwall[icat][isample]->SetFillColor(colerz[icat]);
     stack[istack]->Add(histos[ivol]->hwall[icat][isample]);
   }
  stack[istack]->SetTitle(0);
  stack[istack]->Draw("h");
  stack[istack]->GetXaxis()->SetTitle("Wall [cm]");
  stack[istack]->Draw("h");
  //
  canvas->cd(1);
  istack = 1;
  for (int icat=NCAT-1; icat>=0; icat--){
     histos[ivol]->herec[icat][isample]->SetLineColor(colerz[icat]);
     histos[ivol]->herec[icat][isample]->SetFillColor(colerz[icat]);
     stack[istack]->Add(histos[ivol]->herec[icat][isample]);
  }
  stack[istack]->SetTitle(0);
  stack[istack]->Draw("h");
  stack[istack]->GetXaxis()->SetTitle("E_{rec} [MeV]");
  stack[istack]->Draw("h");

  //
  canvas->cd(5);
  istack = 2;
  for (int icat=NCAT-1; icat>=0; icat--){
     histos[ivol]->htowall[icat][isample]->SetLineColor(colerz[icat]);
     histos[ivol]->htowall[icat][isample]->SetFillColor(colerz[icat]);
     stack[istack]->Add(histos[ivol]->htowall[icat][isample]);
  }
  stack[istack]->SetTitle(0);
  stack[istack]->Draw("h");
  stack[istack]->GetXaxis()->SetTitle("Towall [cm]");
  stack[istack]->Draw("h");
  //
  canvas->cd(2);
  istack = 3;
  for (int icat=NCAT-1; icat>=0; icat--){
     histos[ivol]->hetrue[icat][isample]->SetLineColor(colerz[icat]);
     histos[ivol]->hetrue[icat][isample]->SetFillColor(colerz[icat]);
     stack[istack]->Add(histos[ivol]->hetrue[icat][isample]);
  }
  stack[istack]->SetTitle(0);
  stack[istack]->Draw("h");
  stack[istack]->GetXaxis()->SetTitle("E_{true} [MeV]");
  stack[istack]->Draw("h");
  //
  canvas->cd(4);
  istack = 4;
  for (int icat=NCAT-1; icat>=0; icat--){
     histos[ivol]->hwallv[icat][isample]->SetLineColor(colerz[icat]);
     histos[ivol]->hwallv[icat][isample]->SetFillColor(colerz[icat]);
     stack[istack]->Add(histos[ivol]->hwallv[icat][isample]);
  }
  stack[istack]->SetTitle(0);
  stack[istack]->Draw("h");
  stack[istack]->GetXaxis()->SetTitle("True Wall [cm]");
  // 
  stack[istack]->Draw("h");
  canvas->cd(6);
  istack = 5;
  for (int icat=NCAT-1; icat>=0; icat--){
     histos[ivol]->htowallv[icat][isample]->SetLineColor(colerz[icat]);
     histos[ivol]->htowallv[icat][isample]->SetFillColor(colerz[icat]);
     stack[istack]->Add(histos[ivol]->htowallv[icat][isample]);
  }
  stack[istack]->SetTitle(0);
  stack[istack]->Draw("h");
  stack[istack]->GetXaxis()->SetTitle("True Towall [cm]");
  stack[istack]->Draw("h");


}
*/




//////////////////////////////////////////////////////////////////////////////////////////
void T2KFVStudies::drawPassFail(const char* var,int isamp){

  // communicat
  cout<<"T2KFVStudies: drawing "<<var<<" plot! "<<endl; 
  // canvs set up
  canvas = new TCanvas("cc","cc",900,900);
  canvas->Divide(2,2);

  // draw pass 
  canvas->cd(1);
  histosCompare[0]->drawStack(var,isamp);
  canvas->cd(3);
  histosCompare[1]->drawStack(var,isamp);

  // draw fail 
  canvas->cd(2);
  histosCompare[2]->drawStack(var,isamp);
  canvas->cd(4);
  histosCompare[3]->drawStack(var,isamp);

  // draw with same axes
  double max = histosCompare[0]->stack[isamp]->GetMaximum();
  for (int i=1; i<4; i++){
    double testmax = histosCompare[i]->stack[isamp]->GetMaximum();
    if (testmax > max){
      max = testmax;
    }
  }
  cout<<"T2KFVStudies: Max stack value: "<<max<<endl;
  for (int i=0; i<4; i++){
    histosCompare[i]->stack[isamp]->SetMinimum(0);
    histosCompare[i]->stack[isamp]->SetMaximum(max);
  }

  //
  canvas->cd(1);
  histosCompare[0]->stack[isamp]->Draw("histo");
  //
  canvas->cd(2);
  histosCompare[1]->stack[isamp]->Draw("histo");
  //
  canvas->cd(3);
  histosCompare[2]->stack[isamp]->Draw("histo");
  //
  canvas->cd(4);
  histosCompare[3]->stack[isamp]->Draw("histo");

  // return;
  return;
}



//////////////////////////////////////////////////////////////////////////////////////////
void T2KFVStudies::drawCompare(const char* var){

  // canvs set up
  canvas = new TCanvas("cc","cc",900,900);
  canvas->Divide(2,2);

  // draw first selection
  canvas->cd(1);
  histosCompare[0]->drawStack(var,0);
  canvas->cd(3);
  histosCompare[0]->drawStack(var,1);

  // draw 2nd selection
  canvas->cd(2);
  histosCompare[1]->drawStack(var,0);
  canvas->cd(4);
  histosCompare[1]->drawStack(var,1);

  // draw with same axes
  makeAxesSame(histosCompare[0]->stack[0], histosCompare[1]->stack[0]); 
  makeAxesSame(histosCompare[0]->stack[1], histosCompare[1]->stack[1]); 

  //
  canvas->cd(1);
  histosCompare[0]->stack[0]->Draw("h");
  //
  canvas->cd(2);
  histosCompare[1]->stack[0]->Draw("h");
  //
  canvas->cd(3);
  histosCompare[0]->stack[1]->Draw("h");
  //
  canvas->cd(4);
  histosCompare[1]->stack[1]->Draw("h");

  // return;
  return;
}




//////////////////////////////////////////////////////////////////////////////////////////
void T2KFVStudies::fill(FVStudyHistos* hh){

  // cut-dependant values
  int nu_type; //< expected neutrino flavor for this event
  int icat; //< what catagory (CCQE, CC, CCMisID, NC, Entering)
  int isamp; //< which "sample": 0-> nue CCQE, 1->numu CCQE, 2->nue CC1RPi
  int ifqindex; //< which index to use for fiTQun arrays (1->e, 2->mu)
  float nu_energy;
  //
  if (passedECuts){
    isamp = 0;
    nu_type = 12;
    ifqindex = 1;
    icat = getCatagory(nu_type);
    nu_energy = fqpars.fqenue;
  }
  if (passedMuCuts){
    nu_type = 14;
    isamp = 1;
    ifqindex = 2;
    icat = getCatagory(nu_type);
    nu_energy = fqpars.fqenumu;
  }
  if (passedE1RPiCuts){
    nu_type = 12;
    isamp = 2;
    ifqindex = 1;
    icat = getCatagory(nu_type);
    TVector3 dirrc;
    TVector3 dirtrue;
    dirrc.SetXYZ(mcEvent->fq1rdir[0][1][0],
                 mcEvent->fq1rdir[0][1][1],
                 mcEvent->fq1rdir[0][1][2]);
    nu_energy = calcResPiNuEnergy(mcEvent->fq1rmom[0][1],&dirrc);
    if (nu_energy > 1250.) nu_energy = -1.;
  }


  // fill the histograms
  hh->herec[icat][isamp]->Fill(nu_energy, NormPOT*mcEvent->evtweight);
  hh->herectot[isamp]->Fill(nu_energy, NormPOT*mcEvent->evtweight);
  // fill ccqe residual
  if (TMath::Abs(mcEvent->mode)==1){

    // Erec residual
    float momres =  nu_energy - mcEvent->pmomv[0];
    hh->heres[isamp]->Fill(momres,NormPOT*mcEvent->evtweight);

    // direction residual
    TVector3 dirrc;
    TVector3 dirtrue;
    dirrc.SetXYZ(mcEvent->fq1rdir[0][ifqindex][0],
                 mcEvent->fq1rdir[0][ifqindex][1],
                 mcEvent->fq1rdir[0][ifqindex][2]);
    dirtrue.SetXYZ(mcEvent->dirv[2][0],
                   mcEvent->dirv[2][1],
                   mcEvent->dirv[2][2]);
    float dirres = dirtrue.Angle(dirrc);

    // fill histos
    hh->hdirres[isamp]->Fill(dirres,NormPOT*mcEvent->evtweight);
    hh->heresvwall[isamp]->Fill(mcEvent->fqwall,momres,NormPOT*mcEvent->evtweight);
    hh->heresvtowall[isamp]->Fill(mcEvent->fqtowall,momres,NormPOT*mcEvent->evtweight);
    hh->hfveres[isamp]->Fill(mcEvent->fqtowall,mcEvent->fqwall,momres*NormPOT*mcEvent->evtweight);
  } //< end CCQE block
 
  // fv stuff
  hh->hfvtot[isamp]->Fill(mcEvent->fqtowall,mcEvent->fqwall,NormPOT*mcEvent->evtweight);
  hh->hfvtottrue[isamp]->Fill(mcEvent->towallv[2],mcEvent->wallv,NormPOT*mcEvent->evtweight);
  hh->hwall[icat][isamp]->Fill(mcEvent->fqwall, NormPOT*mcEvent->evtweight);
  hh->hwallv[icat][isamp]->Fill(mcEvent->wallv, NormPOT*mcEvent->evtweight);
  hh->hetrue[icat][isamp]->Fill(mcEvent->pmomv[0], NormPOT*mcEvent->evtweight);
  hh->htowall[icat][isamp]->Fill(mcEvent->fqtowall, NormPOT*mcEvent->evtweight);
  hh->htowallv[icat][isamp]->Fill(mcEvent->towallv[2], NormPOT*mcEvent->evtweight);
  hh->hfv[icat][isamp]->Fill(mcEvent->fqtowall,mcEvent->fqwall,NormPOT*mcEvent->evtweight);
  hh->hpos[icat][isamp]->Fill(TMath::Sqrt(mcEvent->posv[0]*mcEvent->posv[0] + mcEvent->posv[1]*mcEvent->posv[1]),
                                  TMath::Abs(mcEvent->posv[2]),NormPOT*mcEvent->evtweight); 
  // count it
  hh->countEvent(icat, isamp, NormPOT*mcEvent->evtweight); 
  
  //
  return;
}



/*  
//////////////////////////////////////////////////////////////////////////////////////////
void T2KFVStudies::fill(FVStudyHistos* hh){
  
    // get catagory//////////////////////////////
    int icat = 0;
    int isamp = 0;
    TVector3 dirrc;
    TVector3 dirtrue;

    if (passedECuts){
      icat = getCatagory(12);
      isamp = 0;
      hh->herec[icat][isamp]->Fill(fqpars.fqenue, NormPOT*mcEvent->evtweight);
      hh->herectot[0]->Fill(fqpars.fqenue, NormPOT*mcEvent->evtweight);

      // fill ccqe residual
      if (TMath::Abs(mcEvent->mode)==1){
        float momres = (fqpars.fqenue - mcEvent->pmomv[0] );
        hh->heres[0]->Fill(momres,NormPOT*mcEvent->evtweight);
        // direction residual
        dirrc.SetXYZ(mcEvent->fq1rdir[0][1][0],
                     mcEvent->fq1rdir[0][1][1],
                     mcEvent->fq1rdir[0][1][2]);
        dirtrue.SetXYZ(mcEvent->dirv[2][0],
                       mcEvent->dirv[2][1],
                       mcEvent->dirv[2][2]);
        float dirres = dirtrue.Angle(dirrc);
        hh->hdirres[0]->Fill(dirres,NormPOT*mcEvent->evtweight);
        hh->heresvwall[0]->Fill(mcEvent->fqwall,momres,NormPOT*mcEvent->evtweight);
        hh->heresvtowall[0]->Fill(mcEvent->fqtowall,momres,NormPOT*mcEvent->evtweight);
        hh->hfveres[0]->Fill(mcEvent->fqtowall,mcEvent->fqwall,momres*NormPOT*mcEvent->evtweight);
      } //< end CCQE block
 
      // fv stuff
      hh->hfvtot[0]->Fill(mcEvent->fqtowall,mcEvent->fqwall,NormPOT*mcEvent->evtweight);
      hh->hfvtottrue[0]->Fill(mcEvent->towallv[2],mcEvent->wallv,NormPOT*mcEvent->evtweight);
      hh->hwall[icat][isamp]->Fill(mcEvent->fqwall, NormPOT*mcEvent->evtweight);
      hh->hwallv[icat][isamp]->Fill(mcEvent->wallv, NormPOT*mcEvent->evtweight);
      hh->hetrue[icat][isamp]->Fill(mcEvent->pmomv[0], NormPOT*mcEvent->evtweight);
      hh->htowall[icat][isamp]->Fill(mcEvent->fqtowall, NormPOT*mcEvent->evtweight);
      hh->htowallv[icat][isamp]->Fill(mcEvent->towallv[2], NormPOT*mcEvent->evtweight);
      hh->hfv[icat][isamp]->Fill(mcEvent->fqtowall,mcEvent->fqwall,NormPOT*mcEvent->evtweight);
      hh->hpos[icat][isamp]->Fill(TMath::Sqrt(mcEvent->posv[0]*mcEvent->posv[0] + mcEvent->posv[1]*mcEvent->posv[1]),
                                  TMath::Abs(mcEvent->posv[2]),NormPOT*mcEvent->evtweight); 

      // count it
      hh->countEvent(icat, isamp, NormPOT*mcEvent->evtweight); 

    }

    if (passedE1RPiCuts){

      icat = getCatagory(12);
      isamp = 2;
      hh->herec[icat][isamp]->Fill(fqpars.fqenue, NormPOT*mcEvent->evtweight);
      hh->herectot[isamp]->Fill(fqpars.fqenue, NormPOT*mcEvent->evtweight);

      // fill ccqe residual
      if (TMath::Abs(mcEvent->mode)==1){
        float momres = (fqpars.fqenue - mcEvent->pmomv[0] );
        hh->heres[isamp]->Fill(momres,NormPOT*mcEvent->evtweight);
        // direction residual
        dirrc.SetXYZ(mcEvent->fq1rdir[0][1][0],
                     mcEvent->fq1rdir[0][1][1],
                     mcEvent->fq1rdir[0][1][2]);
        dirtrue.SetXYZ(mcEvent->dirv[2][0],
                       mcEvent->dirv[2][1],
                       mcEvent->dirv[2][2]);
        float dirres = dirtrue.Angle(dirrc);
        hh->hdirres[isamp]->Fill(dirres,NormPOT*mcEvent->evtweight);
        hh->heresvwall[isamp]->Fill(mcEvent->fqwall,momres,NormPOT*mcEvent->evtweight);
        hh->heresvtowall[isamp]->Fill(mcEvent->fqtowall,momres,NormPOT*mcEvent->evtweight);
        hh->hfveres[isamp]->Fill(mcEvent->fqtowall,mcEvent->fqwall,momres*NormPOT*mcEvent->evtweight);
      } //< end CCQE block
 
      // fv stuff
      hh->hfvtot[isamp]->Fill(mcEvent->fqtowall,mcEvent->fqwall,NormPOT*mcEvent->evtweight);
      hh->hfvtottrue[isamp]->Fill(mcEvent->towallv[2],mcEvent->wallv,NormPOT*mcEvent->evtweight);
      hh->hwall[icat][isamp]->Fill(mcEvent->fqwall, NormPOT*mcEvent->evtweight);
      hh->hwallv[icat][isamp]->Fill(mcEvent->wallv, NormPOT*mcEvent->evtweight);
      hh->hetrue[icat][isamp]->Fill(mcEvent->pmomv[0], NormPOT*mcEvent->evtweight);
      hh->htowall[icat][isamp]->Fill(mcEvent->fqtowall, NormPOT*mcEvent->evtweight);
      hh->htowallv[icat][isamp]->Fill(mcEvent->towallv[2], NormPOT*mcEvent->evtweight);
      hh->hfv[icat][isamp]->Fill(mcEvent->fqtowall,mcEvent->fqwall,NormPOT*mcEvent->evtweight);
      hh->hpos[icat][isamp]->Fill(TMath::Sqrt(mcEvent->posv[0]*mcEvent->posv[0] + mcEvent->posv[1]*mcEvent->posv[1]),
                                  TMath::Abs(mcEvent->posv[2]),NormPOT*mcEvent->evtweight); 

      // count it
      hh->countEvent(icat, isamp, NormPOT*mcEvent->evtweight); 

    }


    //
    if (passedMuCuts){
      icat = getCatagory(14);
      isamp = 1;
      hh->herec[icat][isamp]->Fill(fqpars.fqenumu,NormPOT*mcEvent->evtweight);
      hh->herectot[1]->Fill(fqpars.fqenumu, NormPOT*mcEvent->evtweight);

      // fill ccqe residual
      if (TMath::Abs(mcEvent->mode)==1){
      float momres = (fqpars.fqenumu - mcEvent->pmomv[0] );
        hh->heres[1]->Fill(momres,NormPOT*mcEvent->evtweight);
        // direction residual
        dirrc.SetXYZ(mcEvent->fq1rdir[0][2][0],
                     mcEvent->fq1rdir[0][2][1],
                     mcEvent->fq1rdir[0][2][2]);
        dirtrue.SetXYZ(mcEvent->dirv[2][0],
                       mcEvent->dirv[2][1],
                       mcEvent->dirv[2][2]);
        float dirres = dirtrue.Angle(dirrc);
        hh->hdirres[1]->Fill(dirres,NormPOT*mcEvent->evtweight);
        hh->hfvtot[1]->Fill(mcEvent->fqtowall,mcEvent->fqwall,NormPOT*mcEvent->evtweight);
        hh->hfvtot[1]->Fill(mcEvent->wallv,mcEvent->towallv[2],NormPOT*mcEvent->evtweight);
        hh->heresvwall[1]->Fill(mcEvent->fqwall,momres,NormPOT*mcEvent->evtweight);
        hh->heresvtowall[1]->Fill(mcEvent->fqtowall,momres,NormPOT*mcEvent->evtweight);
        hh->hfveres[1]->Fill(mcEvent->fqtowall,mcEvent->fqwall,momres*NormPOT*mcEvent->evtweight);
      } //< end CCQE block

     // fv stuff
     hh->hwall[icat][isamp]->Fill(mcEvent->fqwall, NormPOT*mcEvent->evtweight);
     hh->hwallv[icat][isamp]->Fill(mcEvent->wallv, NormPOT*mcEvent->evtweight);
     hh->hfvtottrue[1]->Fill(mcEvent->towallv[2],mcEvent->wallv,NormPOT*mcEvent->evtweight);
     hh->hetrue[icat][isamp]->Fill(mcEvent->pmomv[0], NormPOT*mcEvent->evtweight);
     hh->htowall[icat][isamp]->Fill(mcEvent->fqtowall, NormPOT*mcEvent->evtweight);
     hh->htowallv[icat][isamp]->Fill(mcEvent->towallv[2], NormPOT*mcEvent->evtweight);
     hh->hfv[icat][isamp]->Fill(mcEvent->fqtowall,mcEvent->fqwall,NormPOT*mcEvent->evtweight);
     hh->hpos[icat][isamp]->Fill(TMath::Sqrt(mcEvent->posv[0]*mcEvent->posv[0] + mcEvent->posv[1]*mcEvent->posv[1]),
                                 TMath::Abs(mcEvent->posv[2]),NormPOT*mcEvent->evtweight); 
 
     // count it
     hh->countEvent(icat, isamp, NormPOT*mcEvent->evtweight); 
    }

   
   //
   return;

}
*/


//////////////////////////////////////////////////////////////////////////////////////////
double T2KFVStudies::stackIntegral(THStack* hstack){
  
  double gral = 0.;
  TObjArray* obarr = hstack->GetStack();
  cout<<"adding over "<<obarr->GetEntries()<<" histos "<<endl;
  for (int ih = 0; ih<obarr->GetEntries(); ih++){
    TH1D* hh = (TH1D*)obarr->At(ih);
    cout<<"adding: "<<hh->Integral()<<endl;
    gral += hh->Integral();
  }
 
  return gral;

}



//////////////////////////////////////////////////////////////////////////////////////////
void T2KFVStudies::drawFVSummary(int icat, int ivol, int isamp){

  canvas = new TCanvas("cc","cc",900,900);
  canvas->Divide(2,2);

  canvas->cd(1);
  histos[ivol]->hwallv[icat][isamp]->Draw("h");
  //
  canvas->cd(2);
  histos[ivol]->hwall[icat][isamp]->Draw("h");
  //
  canvas->cd(3);
  histos[ivol]->htowallv[icat][isamp]->Draw("h");
  //
  canvas->cd(4);
  histos[ivol]->htowall[icat][isamp]->Draw("h");

   return;
}


//////////////////////////////////////////////////////////////////////////////////////////
void T2KFVStudies::fillHistos(int ivol){

   fill(histos[ivol]);

   //
   return;
}


 
//////////////////////////////////////////////////////////////////////////////////////////
int T2KFVStudies::passT2KCuts(){
 
   passedCuts = 0;
   passedECuts = 0;
   passedE1RPiCuts= 0;
   passedMuCuts = 0;

  // cuts ////////////////////
  if (mcEvent->fqwall < 0) return 0;
  fqpars.nhitac = mcEvent->nhitac;
  fqpars.fqpid = mcEvent->fq1rnll[0][2]-mcEvent->fq1rnll[0][1];
  fqpars.fqpi0par = mcEvent->fqpi0par;
  fqpars.fqpippar = mcEvent->attribute[2];
  fqpars.fqrcpar = mcEvent->fqrcpar;
  fqpars.fqmommu =mcEvent->fq1rmom[0][2];
  fqpars.fqmome = mcEvent->fq1rmom[0][1];
  fqpars.fqenue = mcEvent->fq1renu[0];
  fqpars.fqenumu = mcEvent->fq1renu[1];
  fqpars.fqnsubev = mcEvent->fqnse;
  fqpars.fqnring = mcEvent->fqmrnring[0];

  if (mcEvent->passecut){
    passedCuts = 1;
    passedECuts = 1;

  }
  if (mcEvent->passmucut){
    passedCuts = 1;
    passedMuCuts = 1;
  }
  if (mcEvent->passe1rpicut){
    passedCuts = 1;
    passedE1RPiCuts = 1;
  }
  

  return passedCuts;
}



//////////////////////////////////////////////////////////////////////////////////////////
int T2KFVStudies::passFVCuts(int fvcode){
 
   // FV cut  
   // code for fv cuts
   // 0 -> no wall or towal cuts
   // 1 -> nominal wall cuts
   // 2 -> inverse nominal cuts
   // 3 -> small wall cut
   // 4 -> entering bg cuts
    if (fvcode==0) {return 1.;}
    else if (fvcode==1){
      if (mcEvent->fqwall > 200.) return 1 ;
    }
    else if (fvcode==2){
      if (mcEvent->fqwall < 200.) return 1;
    }
    else if (fvcode==3){
      if (mcEvent->fqwall < 50.) return 1;
    }
    else if (fvcode==4){
      if (mcEvent->wallv < 0.) return 1;
    }
    else if (fvcode==5){
      if ( (mcEvent->fqwall > 50.) && (mcEvent->fqtowall > 150.)) return 1;
    }
  


  return 0;
}



///////////////////////////////////////////////////////////////////////////////
void T2KFVStudies::compareCuts(double towall1, double wall1, double towall2, double wall2){

  // clear histos
  histosCompare[0]->clearHistos();
  histosCompare[1]->clearHistos();
  histosCompare[2]->clearHistos();
  histosCompare[3]->clearHistos();
  histosCompare[4]->clearHistos();

  int passedcuts1;
  int passedcuts2;

  // get entries
  int N = chMC->GetEntries();

  // loop over evets
//  for (int i=0; i<chMC->GetEntries(); i++){
  for (int i=0; i<Nmax; i++){
//  for (int i=0; i<50000; i++){
    
    // communicate
    if ((i%50000)==0) cout<<i<<endl;

    // read in entry
    chMC->GetEntry(i);

    // does it pass any T2K cuts?
    passedCuts = passT2KCuts();
    if (!passedCuts) continue;

    // cuts for volume 1
    if ( (towall1 < mcEvent->fqtowall) && (wall1 < mcEvent->fqwall) ){
      // fill passed 
      passedcuts1 = 1;
      fill(histosCompare[0]);
    }
    else{
      // fill failed 
      passedcuts1 = 0;
      fill(histosCompare[2]);      
    }

    // cuts fof volume 2
    if ( (towall2 < mcEvent->fqtowall) && (wall2 < mcEvent->fqwall) ){
      fill(histosCompare[1]);
      passedcuts2 = 1;
    }
    else{
      // fill failed 
      fill(histosCompare[3]);      
      passedcuts2 = 0;
    }


    // difference sample 
//    cout<<"diff: "<<(passedcuts1 + passedcuts2)<<endl;
    if ( (passedcuts1 + passedcuts2) == 1 ){
//      cout<<"fillin diff: "<<(passedcuts1 + passedcuts2)<<endl;
      fill(histosCompare[4]);
    }

  } //< end event loop

  // normalize residuals
  histosCompare[0]->normalizeResiduals();
  histosCompare[1]->normalizeResiduals();
  histosCompare[2]->normalizeResiduals();
  histosCompare[3]->normalizeResiduals();
  histosCompare[4]->normalizeResiduals();

  cout<<"#############################"<<endl;
  cout<<"## NU  E ####################"<<endl;
  cout<<"#############################"<<endl;
  cout<<"## PASSED: ("<<towall1<<","<<wall1<<")"<<endl;
  histosCompare[0]->printSummary(0);
  cout<<"## FAILED: ("<<towall1<<","<<wall1<<")"<<endl;
  histosCompare[2]->printSummary(0);
  cout<<"## PASSED: ("<<towall2<<","<<wall2<<")"<<endl;
  histosCompare[1]->printSummary(0);
  cout<<"## FAILED: ("<<towall2<<","<<wall2<<")"<<endl;
  histosCompare[3]->printSummary(0);
  cout<<"#############################"<<endl;
  cout<<"## NU  MU ####################"<<endl;
  cout<<"#############################"<<endl;
  cout<<"## PASSED: ("<<towall1<<","<<wall1<<")"<<endl;
  histosCompare[0]->printSummary(1);
  cout<<"## FAILED: ("<<towall1<<","<<wall1<<")"<<endl;
  histosCompare[2]->printSummary(1);
  cout<<"## PASSED: ("<<towall2<<","<<wall2<<")"<<endl;
  histosCompare[1]->printSummary(1);
  cout<<"## FAILED: ("<<towall2<<","<<wall2<<")"<<endl;
  histosCompare[3]->printSummary(1);
  cout<<"#############################"<<endl;
  cout<<"## NU  E CC1R Pi ############"<<endl;
  cout<<"#############################"<<endl;
  cout<<"## PASSED: ("<<towall1<<","<<wall1<<")"<<endl;
  histosCompare[0]->printSummary(2);
  cout<<"## FAILED: ("<<towall1<<","<<wall1<<")"<<endl;
  histosCompare[2]->printSummary(2);
  cout<<"## PASSED: ("<<towall2<<","<<wall2<<")"<<endl;
  histosCompare[1]->printSummary(2);
  cout<<"## FAILED: ("<<towall2<<","<<wall2<<")"<<endl;
  histosCompare[3]->printSummary(2);


  cout<<"#############################"<<endl;
  cout<<"## DIFF #####################"<<endl;
  cout<<"#############################"<<endl;
  cout<<"- NUMU: "<<endl;
  cout<<"CCQE:      "<<histosCompare[0]->eventCount[0][1]<<" -> "
                    <<histosCompare[1]->eventCount[0][1]<<" ("
                    <<histosCompare[1]->eventCount[0][1]-histosCompare[0]->eventCount[0][1]
                    <<")"<<endl;
  cout<<"CCnQE:     "<<histosCompare[0]->eventCount[1][1]<<" -> "
                    <<histosCompare[1]->eventCount[1][1]<<" ("
                    <<histosCompare[1]->eventCount[1][1]-histosCompare[0]->eventCount[1][1]
                    <<")"<<endl;
  cout<<"CCMisID:   "<<histosCompare[0]->eventCount[2][1]<<" -> "
                    <<histosCompare[1]->eventCount[2][1]<<" ("
                    <<histosCompare[1]->eventCount[2][1]-histosCompare[0]->eventCount[2][1]
                    <<")"<<endl;
  cout<<"NC:        "<<histosCompare[0]->eventCount[3][1]<<" -> "
                    <<histosCompare[1]->eventCount[3][1]<<" ("
                    <<histosCompare[1]->eventCount[3][1]-histosCompare[0]->eventCount[3][1]
                    <<")"<<endl;                    
  cout<<"Entering:  "<<histosCompare[0]->eventCount[4][1]<<" -> "
                    <<histosCompare[1]->eventCount[4][1]<<" ("
                    <<histosCompare[1]->eventCount[4][1]-histosCompare[0]->eventCount[4][1]
                    <<")"<<endl;                    
  cout<<"- NUE: "<<endl;
  cout<<"CCQE:      "<<histosCompare[0]->eventCount[0][0]<<" -> "
                    <<histosCompare[1]->eventCount[0][0]<<" ("
                    <<histosCompare[1]->eventCount[0][0]-histosCompare[0]->eventCount[0][0]
                    <<")"<<endl;
  cout<<"CCnQE:     "<<histosCompare[0]->eventCount[1][0]<<" -> "
                    <<histosCompare[1]->eventCount[1][0]<<" ("
                    <<histosCompare[1]->eventCount[1][0]-histosCompare[0]->eventCount[1][0]
                    <<")"<<endl;
  cout<<"CCMisID:   "<<histosCompare[0]->eventCount[2][0]<<" -> "
                    <<histosCompare[1]->eventCount[2][0]<<" ("
                    <<histosCompare[1]->eventCount[2][0]-histosCompare[0]->eventCount[2][0]
                    <<")"<<endl;
  cout<<"NC:        "<<histosCompare[0]->eventCount[3][0]<<" -> "
                    <<histosCompare[1]->eventCount[3][0]<<" ("
                    <<histosCompare[1]->eventCount[3][0]-histosCompare[0]->eventCount[3][0]
                    <<")"<<endl;                    
  cout<<"Entering:  "<<histosCompare[0]->eventCount[4][0]<<" -> "
                    <<histosCompare[1]->eventCount[4][0]<<" ("
                    <<histosCompare[1]->eventCount[4][0]-histosCompare[0]->eventCount[4][0]
                    <<")"<<endl;                    

  cout<<"- NUE 1R Pi: "<<endl;
  cout<<"CCQE:      "<<histosCompare[0]->eventCount[0][2]<<" -> "
                    <<histosCompare[1]->eventCount[0][2]<<" ("
                    <<histosCompare[1]->eventCount[0][2]-histosCompare[0]->eventCount[0][2]
                    <<")"<<endl;
  cout<<"CCnQE:     "<<histosCompare[0]->eventCount[1][2]<<" -> "
                    <<histosCompare[1]->eventCount[1][2]<<" ("
                    <<histosCompare[1]->eventCount[1][2]-histosCompare[0]->eventCount[1][2]
                    <<")"<<endl;
  cout<<"CCMisID:   "<<histosCompare[0]->eventCount[2][2]<<" -> "
                    <<histosCompare[1]->eventCount[2][2]<<" ("
                    <<histosCompare[1]->eventCount[2][2]-histosCompare[0]->eventCount[2][2]
                    <<")"<<endl;
  cout<<"NC:        "<<histosCompare[0]->eventCount[3][2]<<" -> "
                    <<histosCompare[1]->eventCount[3][2]<<" ("
                    <<histosCompare[1]->eventCount[3][2]-histosCompare[0]->eventCount[3][2]
                    <<")"<<endl;                    
  cout<<"Entering:  "<<histosCompare[0]->eventCount[4][2]<<" -> "
                    <<histosCompare[1]->eventCount[4][2]<<" ("
                    <<histosCompare[1]->eventCount[4][2]-histosCompare[0]->eventCount[4][2]
                    <<")"<<endl;                    

  //
  return;
}



///////////////////////////////////////////////////////////////////////////////
void T2KFVStudies::drawSummary(int ivol, int isample){


  // setup stacks
  for (int is=0; is<10; is++){
    stack[is] = new THStack(Form("hs%d",is),Form("hstack%d",is));
  }

  // color scheme
  int colerz[NCAT] = {4,7,2,1,6};
 
  // canvas
  canvas = new TCanvas("cc","cc",800,900);
  canvas->Divide(2,3);

  // drawing
  canvas->cd(3);
  int istack = 0;
  for (int icat=NCAT-1; icat>=0; icat--){
     histos[ivol]->hwall[icat][isample]->SetLineColor(colerz[icat]);
     histos[ivol]->hwall[icat][isample]->SetFillColor(colerz[icat]);
     stack[istack]->Add(histos[ivol]->hwall[icat][isample]);
   }
  stack[istack]->SetTitle(0);
  stack[istack]->Draw("h");
  stack[istack]->GetXaxis()->SetTitle("Wall [cm]");
  stack[istack]->Draw("h");
  //
  canvas->cd(1);
  istack = 1;
  for (int icat=NCAT-1; icat>=0; icat--){
     histos[ivol]->herec[icat][isample]->SetLineColor(colerz[icat]);
     histos[ivol]->herec[icat][isample]->SetFillColor(colerz[icat]);
     stack[istack]->Add(histos[ivol]->herec[icat][isample]);
  }
  stack[istack]->SetTitle(0);
  stack[istack]->Draw("h");
  stack[istack]->GetXaxis()->SetTitle("E_{rec} [MeV]");
  stack[istack]->Draw("h");

  //
  canvas->cd(5);
  istack = 2;
  for (int icat=NCAT-1; icat>=0; icat--){
     histos[ivol]->htowall[icat][isample]->SetLineColor(colerz[icat]);
     histos[ivol]->htowall[icat][isample]->SetFillColor(colerz[icat]);
     stack[istack]->Add(histos[ivol]->htowall[icat][isample]);
  }
  stack[istack]->SetTitle(0);
  stack[istack]->Draw("h");
  stack[istack]->GetXaxis()->SetTitle("Towall [cm]");
  stack[istack]->Draw("h");
  //
  canvas->cd(2);
  istack = 3;
  for (int icat=NCAT-1; icat>=0; icat--){
     histos[ivol]->hetrue[icat][isample]->SetLineColor(colerz[icat]);
     histos[ivol]->hetrue[icat][isample]->SetFillColor(colerz[icat]);
     stack[istack]->Add(histos[ivol]->hetrue[icat][isample]);
  }
  stack[istack]->SetTitle(0);
  stack[istack]->Draw("h");
  stack[istack]->GetXaxis()->SetTitle("E_{true} [MeV]");
  stack[istack]->Draw("h");
  //
  canvas->cd(4);
  istack = 4;
  for (int icat=NCAT-1; icat>=0; icat--){
     histos[ivol]->hwallv[icat][isample]->SetLineColor(colerz[icat]);
     histos[ivol]->hwallv[icat][isample]->SetFillColor(colerz[icat]);
     stack[istack]->Add(histos[ivol]->hwallv[icat][isample]);
  }
  stack[istack]->SetTitle(0);
  stack[istack]->Draw("h");
  stack[istack]->GetXaxis()->SetTitle("True Wall [cm]");
  // 
  stack[istack]->Draw("h");
  canvas->cd(6);
  istack = 5;
  for (int icat=NCAT-1; icat>=0; icat--){
     histos[ivol]->htowallv[icat][isample]->SetLineColor(colerz[icat]);
     histos[ivol]->htowallv[icat][isample]->SetFillColor(colerz[icat]);
     stack[istack]->Add(histos[ivol]->htowallv[icat][isample]);
  }
  stack[istack]->SetTitle(0);
  stack[istack]->Draw("h");
  stack[istack]->GetXaxis()->SetTitle("True Towall [cm]");
  stack[istack]->Draw("h");


  return;
}


///////////////////////////////////////////////////////////////////////////////
int T2KFVStudies::getCatagory(int wantnutype){

    // is dead region?
    if (mcEvent->wallv < 0.) { return 4;}

    // is NC?
    if (TMath::Abs(mcEvent->mode)>=30){return 3;}

    // is Mis ID?
    if (TMath::Abs(mcEvent->ipnu[0])!=wantnutype) {return 2;}

    // is CCQE?
    if (TMath::Abs(mcEvent->mode)<10){return 0;}

    // CCnQE 
    if ((TMath::Abs(mcEvent->mode)>=10)&&(TMath::Abs(mcEvent->mode)<30)) {return 1;}
 
    cout<<"!! no catagory found"<<endl;
    return -1;
}
    


///////////////////////////////////////////////////////////////////////////////
// 
void T2KFVStudies::findPoorTowall(){


  // histo setup
//  hpoorfv[0] = new TH2FV("hpoorfv_nue",-1,30,0,500,30,0,500);
//  hpoorfv[1] = new TH2FV("hpoorfv_numu",-1,30,0,500,30,0,500);
  hpoorfv[0] = new TH2D("hpoorfv_nue","hpoorfv_nue",50,0,5000,50,0,500);
  hpoorfv[1] = new TH2D("hpoorfv_numu","hpoorfv_numu",50,0,5000,50,0,500);
  hpoor[0] = new TH1D("hpoorwall_nue","npoorwall_nue",50,0,500);
  hpoor[1] = new TH1D("hpoorwall_numu","npoorwall_nue",50,0,500);
  hpoortw[0] = new TH1D("hpoortowall_nue","npoortowall_nue",50,0,5000);
  hpoortw[1] = new TH1D("hpoortowall_numu","npoortowall_nue",50,0,5000);
  hpoorlmom[0] = new TH1D("hpoorlmom_nue","npoorlmim_nue",50,0,2000);
  hpoorlmom[1] = new TH1D("hpoorlmom_numu","npoorlmom_numu",50,0,2000);

  int N = chMC->GetEntries();
  for (int i=0; i<chMC->GetEntries(); i++){

    if ((i%50000)==0) cout<<i<<endl;

    // read in entry
    chMC->GetEntry(i);

    // define poor towall
    int flgpoor = 0;
    if (TMath::Abs(mcEvent->fqtowall - mcEvent->towallv[2])>200.){
      if (mcEvent->towallv[2]<200. && TMath::Abs(mcEvent->mode)<10){
        flgpoor = 1;
      }
    }

    if (flgpoor){
      // fill
      if (mcEvent->passecut){
         hpoorfv[0]->Fill(mcEvent->fqtowall, mcEvent->fqwall, mcEvent->evtweight);  
         hpoor[0]->Fill(mcEvent->fqwall, mcEvent->evtweight);  
         hpoortw[0]->Fill(mcEvent->fqtowall, mcEvent->evtweight);  
         hpoorlmom[0]->Fill(mcEvent->fq1rmom[0][1], mcEvent->evtweight);  
      }
      else if (mcEvent->passmucut){
         hpoorfv[1]->Fill(mcEvent->fqtowall, mcEvent->fqwall, mcEvent->evtweight);  
         hpoor[1]->Fill(mcEvent->fqwall, mcEvent->evtweight);  
         hpoortw[1]->Fill(mcEvent->fqtowall, mcEvent->evtweight);  
         hpoorlmom[1]->Fill(mcEvent->fq1rmom[0][2], mcEvent->evtweight);  
      }
    }

    // loop over volumes and fill
    // FV cut  
    // code for fv cuts
    // 0 -> no wall or towal cuts
    // 1 -> nominal wall cuts
    // 2 -> inverse nominal cuts
    // 3 -> small wall cut
    // 4 -> entering bg cuts

  }

  return;

}


///////////////////////////////////////////////////////////////////////////////
// not that useful anymore...
void T2KFVStudies::loop(){

  int N = chMC->GetEntries();
  for (int i=0; i<chMC->GetEntries(); i++){
//  for (int i=0; i<50000; i++){

    if ((i%50000)==0) cout<<i<<endl;
    

    // read in entry
    chMC->GetEntry(i);

    // does it pass any T2K cuts?
    passedCuts = passT2KCuts();
    if (!passedCuts) continue;

    // loop over volumes and fill
    // FV cut  
    // code for fv cuts
    // 0 -> no wall or towal cuts
    // 1 -> nominal wall cuts
    // 2 -> inverse nominal cuts
    // 3 -> small wall cut
    // 4 -> entering bg cuts
    for (int ivol=0; ivol<NVOLS; ivol++){
      if (passFVCuts(ivol)){
        fillHistos(ivol);
      }
    }

  }

  return;

}



///////////////////////////////////////////////////////////////////////////////
void T2KFVStudies::init(){

  // branches
  chMC->SetBranchStatus("*",0);
  chMC->SetBranchStatus("*v",1);
  chMC->SetBranchStatus("fq*",1);
  chMC->SetBranchStatus("mode",1);
  chMC->SetBranchStatus("wall*",1);
  chMC->SetBranchStatus("towall*",1);
  chMC->SetBranchStatus("att*",1);
  chMC->SetBranchStatus("ip*",1);
  chMC->SetBranchStatus("ev*",1);
  chMC->SetBranchStatus("nhitac",1);
  chMC->SetBranchStatus("pass*",1);

  // histogram containters
  for (int ivol=0; ivol<NVOLS; ivol++){
    histos[ivol] = new FVStudyHistos(Form("volume_%d",ivol));
  }
  histosCompare[0] = new FVStudyHistos("Passed_cut1"); 
  histosCompare[1] = new FVStudyHistos("Passed_cut2"); 
  histosCompare[2] = new FVStudyHistos("Failed_cut1"); 
  histosCompare[3] = new FVStudyHistos("Failed_cut2"); 
  histosCompare[4] = new FVStudyHistos("Difference"); 


  /*
  // wall
  TH1D* wallseed = new TH1D("wallseed","wallseed",50,0,1000);
  for (int icat=0; icat<NCAT; icat++){
    for (int isamp=0; isamp<NSAMP; isamp++){
      TString hname = Form("hwall_samp%d_cat%d",isamp,icat);
      hwall[icat][isamp] = (TH1D*)wallseed->Clone(hname.Data());
    }
  }

  // erec
  TH1D* erecseed = new TH1D("erec","erec",50,0,2000);
  for (int icat=0; icat<NCAT; icat++){
    for (int isamp=0; isamp<NSAMP; isamp++){
      TString hname = Form("herec_samp%d_cat%d",isamp,icat);
      herec[icat][isamp] = (TH1D*)erecseed->Clone(hname.Data());
    }
  }
  */

//  hWallAll[0] = (TH1D*)wallseed->Clone("hwall_all_e");
//  hWallCCQE[0] = (TH1D*)wallseed->Clone("hwall_ccqe_e");
//  hWallCCnQE[0] = (TH1D*)wallseed->Clone("hwall_ccnqe_e");
//  hWallCCMisID[0] = (TH1D*)wallseed->Clone("hwall_ccmisid_e");
//  hWallNC[0] = (TH1D*)wallseed->Clone("hwall_nc_e");
  // 
//  hWallAll[1] = (TH1D*)wallseed->Clone("hwall_all_mu");
//  hWallCCQE[1] = (TH1D*)wallseed->Clone("hwall_ccqe_mu");
//  hWallCCnQE[1] = (TH1D*)wallseed->Clone("hwall_ccnqe_mu");
//  hWallCCMisID[1] = (TH1D*)wallseed->Clone("hwall_ccmisid_mu");
//  hWallNC[1] = (TH1D*)wallseed->Clone("hwall_nc_mu");


}



///////////////////////////////////////////////////////////////////////////////
T2KFVStudies::T2KFVStudies(TChain* ch){

  chMC = ch;
  Nmax = chMC->GetEntries();
  mcEvent = new fqProcessedEvent(ch);
  NormPOT = 1.0;
  init();

}




#endif
