
#include <iostream>
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TString.h"
#include "TGraph.h"
#include "TFile.h"
#include "TObjArray.h"
#include "Mar13DataReader.C"
#include "Mar13Reader.C"

#define CINTMODE

using namespace std;


class makePlots{
  public:

  makePlots(TChain *chmc, TChain *chdata);

  TChain* chMC;
  TChain* chData;
  TChain* chRFG;
  Mar13Reader* mcEvent;
  Mar13DataReader* dataEvent;
  TH1D* hh[10][2];
  double Norm;
  double getWeight();
  int flgRFG;
  float rfgwgt;

//  void singlefileloop(const char* filename);
//  TString* getRFGFileName(const char* mcfilename);
  void loop();
  void findFVEvents(float wmin, float wmax, float twmin, float twmax);
  void init(); 
  void compare(int ih);

};

/*
void makePlots::getRFGFileName(const char* mcfilename, TChain* rfgWeightChain){

  // name of found file
  TString foundfile = "";

  // name of mc file
  TString mcfile = rootfile;

  // remove ".root" from file name
  mcfile.Remove(mcfile.Sizeof()-6,5);

  // remove all but last 3 chars
  mcfile.Remove(0,mcfile.Sizeof()-4);
  // ^ should just be a number now

  // get list of files
  TObjArray* oblist = rfgWeightChain->GetListOfFiles();

  // look for this number
  for (int ifile=0; ifile<oblist->GetSize(); ifile++){
    TString testname = oblist->At(ifile)->GetTitle();
    if (testname.Contains(mcfile.Data())){
      foundfile = testname.Data();
      cout<<"Found RFG file: "<<testname.Data()<<" for: "<<rootfile<<endl;
    }
  }

  //
  return foundfile;


}
*/

void makePlots::compare(int ih){
 hh[ih][0]->SetLineColor(kRed);
 hh[ih][0]->Draw("h");
 hh[ih][1]->Draw("samee");
}


void makePlots::init(){
 
 hh[0][0] =  new TH1D("hpid","hpid",30,0,5);
 hh[0][1] =  new TH1D("hpid_data","hpid_data",30,0,5);
 //
 hh[1][0] =  new TH1D("hnring","hnring",10,0,10);
 hh[1][1] =  new TH1D("hnring_data","hnring_dat",10,0,10);

 return;
}



double makePlots::getWeight(){

  double ww = 1.0;

  // get oscillation weight
  ww *=  mcEvent->oscwgt;

  // get flux weight
  ww *= ( 0.65*mcEvent->flxh11[0] + 0.35*mcEvent->flxh11[2] )/mcEvent->flxh06[1];

  if (flgRFG){
    ww *= rfgwgt;
  }

  return ww;
}

/*
void makePlots::singlefileloop(const char* fname){

  // get tree
  TFile* mcfile = new TFile(fname);
  TTree* mctree = (TTree*)mcfile->Get("h1");
  // get reader for it
  mcEvent = new Mar13Reader(mctree);  

  // get rfg weight
  TString rfgfile = getRFGFileName(fname);
  TTree* rfgtree = (TTree*)rfgfile->Get("weightstree");
  // get branch
  float rfgwgt;
  rfgtree->SetBranchAddress("fWeight",&rfgwgt);

  // loop over mc
  int NN = chMC->GetEntries();
  for (int iev=0; iev<NN; iev++){

    if ((iev%5000)==0) cout<<iev<<endl;

    chMC->GetEntry(iev);
    double ww = getWeight();

    // multiply by rfgwgt
    ww *= rfgwgt;

    // selections //////////////////
    if (mcEvent->nhitac >= 16) continue;
    if (mcEvent->evis < 30) continue;
    if (mcEvent->wall < 200) continue;
    ////////////////////////////////

    // fill
    hh[1][0]->Fill(mcEvent->nring,ww*Norm);
   
 }

  return;
}
*/



////////////////////////////////////////////////////////////
void makePlots::findFVEvents(float wmin, float wmax, float twmin, float twmax){

// chMC->SetBranchStatus("*",0);
// chMC->SetBranchStatus("nhitac",1);
// chMC->SetBranchStatus("osc*",1);
// chMC->SetBranchStatus("flx*",1);
// chMC->SetBranchStatus("nring",1);
// chMC->SetBranchStatus("wall",1);
// chMC->SetBranchStatus("evis",1);

 // loop over mc
 
 if (flgRFG){
   chRFG->SetBranchAddress("fWeight",&rfgwgt);
 }

 int NN = chMC->GetEntries();
 for (int iev=0; iev<NN; iev++){
   if ((iev%5000)==0) cout<<iev<<endl;
   chMC->GetEntry(iev);
   if (flgRFG){
     chRFG->GetEntry(iev);
   }
   double ww = getWeight();

   // selections //////////////////
   if (mcEvent->nhitac >= 16) continue;
   if (mcEvent->evis < 30) continue;
   if (mcEvent->wall < 200) continue;
   ////////////////////////////////

   // fill
   hh[1][0]->Fill(mcEvent->nring,ww*Norm);


 }

 
 // loop over dat
 for (int iev=0; iev<chData->GetEntries(); iev++){
   if ((iev%5000)==0) cout<<iev<<endl;
   chData->GetEntry(iev);
   // selections //////////////////
   if (dataEvent->nhitac >= 16) continue;
   if (dataEvent->evis < 30) continue;
   if (dataEvent->wall < 200) continue;
   ////////////////////////////////

   // fill
   hh[1][1]->Fill(dataEvent->nring);
   
 }



 return;
}




////////////////////////////////////////////////////////////
void makePlots::loop(){

// chMC->SetBranchStatus("*",0);
// chMC->SetBranchStatus("nhitac",1);
// chMC->SetBranchStatus("osc*",1);
// chMC->SetBranchStatus("flx*",1);
// chMC->SetBranchStatus("nring",1);
// chMC->SetBranchStatus("wall",1);
// chMC->SetBranchStatus("evis",1);

 // loop over mc
 
 if (flgRFG){
   chRFG->SetBranchAddress("fWeight",&rfgwgt);
 }

 int NN = chMC->GetEntries();
 for (int iev=0; iev<NN; iev++){
   if ((iev%5000)==0) cout<<iev<<endl;
   chMC->GetEntry(iev);
   if (flgRFG){
     chRFG->GetEntry(iev);
   }
   double ww = getWeight();

   // selections //////////////////
   if (mcEvent->nhitac >= 16) continue;
   if (mcEvent->evis < 30) continue;
   if (mcEvent->wall < 200) continue;
   ////////////////////////////////

   // fill
   hh[1][0]->Fill(mcEvent->nring,ww*Norm);


 }

 
 // loop over dat
 for (int iev=0; iev<chData->GetEntries(); iev++){
   if ((iev%5000)==0) cout<<iev<<endl;
   chData->GetEntry(iev);
   // selections //////////////////
   if (dataEvent->nhitac >= 16) continue;
   if (dataEvent->evis < 30) continue;
   if (dataEvent->wall < 200) continue;
   ////////////////////////////////

   // fill
   hh[1][1]->Fill(dataEvent->nring);
   
 }



 return;
}



makePlots::makePlots(TChain *chmc, TChain *chdata){
  
  chMC = chmc;
  chData = chdata;
 
  mcEvent = new Mar13Reader(chMC);
  dataEvent= new Mar13DataReader(chData);

  double livetime = 2166.5;
  double fMCyr = (double)chMC->GetNtrees();
  double factor = livetime/365.25/fMCyr;
  Norm = factor;
  flgRFG = 0;
  init();
}




