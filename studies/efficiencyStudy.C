#ifndef EFFICIENCYSTUDY_C
#define EFFICIENCYSTUDY_C

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


///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
class efficiencyStudy{
 
  public:


  efficiencyStudy(TChain *chm, TChain* chd);

  int NCuts;

  TChain* chMC;
  TChain* chDat;

  fqProcessedEvent* mcEvent;
  fqProcessedEvent* datEvent;

  TH2FV* hFVTotal[2];
  TH2FV* hFVCore[2];

  TH2FV* hFVTotalData[2];
  TH2FV* hFVCoreData[2];
 
  TH2FV* hEff[2];
  TH2FV* hEffData[2];

  TH2FV* hRelativeEff[2];
  TH2FV* hRelativeEffData[2];

  TH2FV* hEffCompare[2];
  TH2FV* hRelEffCompare[2];

  TH2FV* hFVUncMap[2];

  int getClass(double pid, int nsubev, double evis);
  int isCore(int iclass, int nring, double pi0par, double pippar);

  void fillHistos();

  void calcEfficiency();

  void initHistos(int bintype=0);
};
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////
void efficiencyStudy::calcEfficiency(){

  int nbinsmax = hEff[0]->GetNumberOfBins();

  // efficiency
  for (int ibin=1; ibin<=nbinsmax; ibin++){
    for (int iclass=0; iclass<2; iclass++){
      double ntot = hFVTotal[iclass]->GetBinContent(ibin);
      double ncore = hFVCore[iclass]->GetBinContent(ibin);
      double eff;
      if (ntot==0) eff=0.;
      else{
        eff = 100.*ncore/ntot;
      }
      hEff[iclass]->SetBinContent(ibin,eff);
      //
      ntot = hFVTotalData[iclass]->GetBinContent(ibin);
      ncore = hFVCoreData[iclass]->GetBinContent(ibin);
      if (ntot==0) eff=0.;
      else{
        eff = 100.*ncore/ntot;
      }
      hEffData[iclass]->SetBinContent(ibin,eff);
    }
  }


  // relativeefficiency
  for (int ibin=1; ibin<=nbinsmax; ibin++){
    for (int iclass=0; iclass<2; iclass++){
      
      double centerval = hEff[iclass]->GetBinContent(nbinsmax);
      double effdiff =  hEff[iclass]->GetBinContent(ibin) - centerval;
      double frac;
      if (centerval==0){
        frac = 0.;
      }
      else{
//        frac = (hEff[iclass]->GetBinContent(ibin)/centerval - 1.);
        frac = 100.*(hEff[iclass]->GetBinContent(ibin)/centerval);
      }
      hRelativeEff[iclass]->SetBinContent(ibin,frac);
      //
      centerval = hEffData[iclass]->GetBinContent(nbinsmax);
      effdiff =  hEffData[iclass]->GetBinContent(ibin) - centerval;
      if (centerval==0){
        frac = 0.;
      }
      else{
//        frac = (hEff[iclass]->GetBinContent(ibin)/centerval - 1.);
        frac = 100.*(hEffData[iclass]->GetBinContent(ibin)/centerval);
      }
      hRelativeEffData[iclass]->SetBinContent(ibin,frac);
    }
  }

  // data/MC compare
  for (int ibin=1; ibin<=nbinsmax; ibin++){
    for (int iclass=0; iclass<2; iclass++){
    
      double mceff = hEff[iclass]->GetBinContent(ibin);
      double dateff = hEffData[iclass]->GetBinContent(ibin);
      double frac;
      if (dateff==0){
        frac =0; 
      }
      else{
        frac = ((mceff/dateff) - 1.)*100.;
      }
      hEffCompare[iclass]->SetBinContent(ibin,frac);
      //
      mceff = hRelativeEff[iclass]->GetBinContent(ibin);
      dateff = hRelativeEffData[iclass]->GetBinContent(ibin);
      if (dateff==0){
        frac =0; 
      }
      else{
        frac = ((mceff/dateff) -1.)*100.;
      }
      hRelEffCompare[iclass]->SetBinContent(ibin,frac);
 
    }

  }

  // uncertainty map
  for (int ibin=1; ibin<=nbinsmax; ibin++){
    for (int iclass=0; iclass<2; iclass++){
      double binc = TMath::Abs(hRelEffCompare[iclass]->GetBinContent(ibin));
      hFVUncMap[iclass]->SetBinContent(ibin, TMath::Abs(binc)/100. ); 
    }
  }


  // set ranges
  for (int iclass=0; iclass<2; iclass++){
    hEff[iclass]->SetMinimum(0.);
    hEff[iclass]->SetMaximum(100.);
    hEffData[iclass]->SetMinimum(0.);
    hEffData[iclass]->SetMaximum(100.);
    hRelativeEff[iclass]->SetMinimum(-20.);
    hRelativeEff[iclass]->SetMaximum(20.);
    hRelativeEffData[iclass]->SetMinimum(-20.);
    hRelativeEffData[iclass]->SetMaximum(20.);
    hEffCompare[iclass]->SetMinimum(-20.);
    hEffCompare[iclass]->SetMaximum(20.);
    hRelEffCompare[iclass]->SetMinimum(-20.);
    hRelEffCompare[iclass]->SetMaximum(20.);
    hFVUncMap[iclass]->SetMinimum(0.);
    hFVUncMap[iclass]->SetMaximum(0.15);
  }

  return;
}


///////////////////////////////////////////////////////////////////
int efficiencyStudy::isCore(int iclass, int nring, double pi0par, double pippar){

  if (NCuts==2){
   if (iclass==0){
     if (nring==1 && pi0par<=0.) return 1;
     else {return 0;}
   }
   if (iclass==1){
     if (nring==1 && pippar<=0.) return 1;
     else {return 0;}
   }
  }
  //
  else if (NCuts==1){
   if (iclass==0){
     if (nring==1) return 1;
     else {return 0;}
   }
   if (iclass==1){
     if (nring==1) return 1;
     else {return 0;}
   }
  }
  //
  if (NCuts==0){
   if (iclass==0){
     if (pi0par<=0.) return 1;
     else {return 0;}
   }
   if (iclass==1){
     if (pippar<=0.) return 1;
     else {return 0;}
   }
  }

  

  //
  return 0;
}



///////////////////////////////////////////////////////////////////
int efficiencyStudy::getClass(double pid, int nsubev, double evis){
 
  // CODE:
  // 0 -> E class
  // 1 -> Mu class
  // -1 -> No class (ha)

  if ( pid>0 && nsubev==1 && evis>100.){
    return 0;
  }
  if ( pid<0 && nsubev<=2 && evis>100.){
    return 1;
  }
  return -1.;
}


///////////////////////////////////////////////////////////////////
void efficiencyStudy::fillHistos(){

  int NevMC = chMC->GetEntries();
  int NevDat = chDat->GetEntries();

  // fill MC
  for (int iev=0; iev<NevMC; iev++){
    
    // read entry
    chMC->GetEntry(iev);
   
    // find class
    int nclass = getClass(mcEvent->attribute[0],(int)mcEvent->fqnse,mcEvent->fq1rmom[0][1]);
    if (nclass<0) continue;
   
    // fill histo 
    hFVTotal[nclass]->Fill(mcEvent->fqtowall, mcEvent->fqwall, mcEvent->evtweight);
   
    // but is it core?
    if (isCore(nclass, mcEvent->fqmrnring[0], mcEvent->attribute[1], mcEvent->attribute[2])){
      hFVCore[nclass]->Fill(mcEvent->fqtowall, mcEvent->fqwall, mcEvent->evtweight);
    }
  }


  // fill data
  for (int iev=0; iev<NevDat; iev++){
    
    // read entry
    chDat->GetEntry(iev);

    // find class
    int nclass = getClass(datEvent->attribute[0],(int)datEvent->fqnse,datEvent->fq1rmom[0][1]);
    if (nclass<0) continue;
   
    // fill histo 
    hFVTotalData[nclass]->Fill(datEvent->fqtowall, datEvent->fqwall);
   
    // but is it core?
    if (isCore(nclass, datEvent->fqmrnring[0], datEvent->attribute[1], datEvent->attribute[2])){
      hFVCoreData[nclass]->Fill(datEvent->fqtowall, datEvent->fqwall);
    }

  }


  return;
}


///////////////////////////////////////////////////////////////////
void efficiencyStudy::initHistos(int bintype){

  // preset binning
  if (bintype>=0){
  for (int i=0; i<2; i++){
    TString hname; 
    //
    hname = Form("total_selection%d",i);
    hFVTotal[i] = new TH2FV(hname.Data(),bintype);
    hname = Form("core_selection%d",i);
    hFVCore[i] = new TH2FV(hname.Data(),bintype);
    hname = Form("total_eff_selection%d",i);
    hEff[i] = new TH2FV(hname.Data(),bintype);
    hname = Form("relative_eff_selection%d",i);
    hRelativeEff[i] = new TH2FV(hname.Data(),bintype);
    //
    hname = Form("total_data_selection%d",i);
    hFVTotalData[i] = new TH2FV(hname.Data(),bintype);
    hname = Form("core_data_selection%d",i);
    hFVCoreData[i] = new TH2FV(hname.Data(),bintype);
    hname = Form("total_eff_selection%d",i);
    hEffData[i] = new TH2FV(hname.Data(),bintype);
    hname = Form("relative_eff_selection%d",i);
    hRelativeEffData[i] = new TH2FV(hname.Data(),bintype);
    // 
    hname = Form("compare_eff_selection%d",i);
    hEffCompare[i] = new TH2FV(hname.Data(),bintype);
    hname = Form("compare_relative_eff_selection%d",i);
    hRelEffCompare[i] = new TH2FV(hname.Data(),bintype);
    hname = "RelativeEfficiencyUncNuE";
    hFVUncMap[0] = new TH2FV(hname.Data(),bintype);
    hname = "RelativeEfficiencyUncNuMu";
    hFVUncMap[1] = new TH2FV(hname.Data(),bintype);
  }
  }

  // other binning
  else if (bintype<0){
  int nbins = 10;
  double xmin = 0;
  double xmax = 2000;
  for (int i=0; i<2; i++){
    TString hname; 
    //
    hname = Form("total_selection%d",i);
    hFVTotal[i] = new TH2FV(hname.Data(),bintype,nbins,xmin,xmax,nbins,xmin,xmax);
    hname = Form("core_selection%d",i);
    hFVCore[i] = new TH2FV(hname.Data(),bintype,nbins,xmin,xmax,nbins,xmin,xmax);
    hname = Form("total_eff_selection%d",i);
    hEff[i] = new TH2FV(hname.Data(),bintype,nbins,xmin,xmax,nbins,xmin,xmax);
    hname = Form("relative_eff_selection%d",i);
    hRelativeEff[i] = new TH2FV(hname.Data(),bintype,nbins,xmin,xmax,nbins,xmin,xmax);
    //
    hname = Form("total_data_selection%d",i);
    hFVTotalData[i] = new TH2FV(hname.Data(),bintype,nbins,xmin,xmax,nbins,xmin,xmax);
    hname = Form("core_data_selection%d",i);
    hFVCoreData[i] = new TH2FV(hname.Data(),bintype,nbins,xmin,xmax,nbins,xmin,xmax);
    hname = Form("total_eff_selection%d",i);
    hEffData[i] = new TH2FV(hname.Data(),bintype,nbins,xmin,xmax,nbins,xmin,xmax);
    hname = Form("relative_eff_selection%d",i);
    hRelativeEffData[i] = new TH2FV(hname.Data(),bintype,nbins,xmin,xmax,nbins,xmin,xmax);
    // 
    hname = Form("compare_eff_selection%d",i);
    hEffCompare[i] = new TH2FV(hname.Data(),bintype,nbins,xmin,xmax,nbins,xmin,xmax);
    hname = Form("compare_relative_eff_selection%d",i);
    hRelEffCompare[i] = new TH2FV(hname.Data(),bintype,nbins,xmin,xmax,nbins,xmin,xmax);
    //
    hname = "RelativeEfficiencyUncNuE";
    hFVUncMap[0] = new TH2FV(hname.Data(),bintype,nbins,xmin,xmax,nbins,xmin,xmax);
    hname = "RelativeEfficiencyUncNuMu";
    hFVUncMap[1] = new TH2FV(hname.Data(),bintype,nbins,xmin,xmax,nbins,xmin,xmax);
  }
  }

  return;

}




///////////////////////////////////////////////////////////////////
efficiencyStudy::efficiencyStudy(TChain *chm, TChain *chd){

  NCuts = 2;
  chMC=chm;
  chDat=chd;
  mcEvent = new fqProcessedEvent(chMC);
  datEvent = new fqProcessedEvent(chDat);

  //
  return;
}

















#endif
