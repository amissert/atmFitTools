
#include <iostream>
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TString.h"
#include "TGraph.h"
#include "TList.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TObjArray.h"
#include "../fqProcessedEvent.cxx"

using namespace std;


class selectData{

  public:

  selectData(TChain* ch);
  TChain* chd;
  fqProcessedEvent* datEvent;

  void selectFV(float w1, float w2, float tw1, float tw2);
  void printEvent(int evtnum);

};


void selectData::printEvent(int evtnum){

  int foundevt = -1;
  for (int iev=0; iev<chd->GetEntries(); iev++){
    chd->GetEntry(iev);
    if (datEvent->nev == evtnum){
      foundevt = iev;
      cout<<"found it! "<<iev<<endl;
      break;
    }
  }

  if (foundevt>0){
    // print info
    chd->GetEntry(foundevt);
    cout<<endl;
    cout<<endl;
    cout<<endl;
    cout<<endl;
    cout<<endl;
    cout<<endl;
    cout<<endl;
    cout<<"---------------------------------------------------"<<endl;
    cout<<"FQWALL:   "<<datEvent->fqwall<<endl;
    cout<<"FQTOWALL: "<<datEvent->fqtowall<<endl;
    cout<<"FQEMOM:   "<<datEvent->fq1rmom[0][1]<<endl;
    cout<<"FQNRING:  "<<datEvent->fqmrnring[0]<<endl;
    cout<<"FQRCPAR:  "<<datEvent->fqrcpar<<endl;
    if (datEvent->fqmrnring[0]>1){
      cout<<"  FQMRMOM1:    "<<datEvent->fqmrmom[datEvent->best2RID][0]<<endl;
      cout<<"  FQMRMOM2:    "<<datEvent->fqmrmom[datEvent->best2RID][1]<<endl;
//      cout<<"  FQMRCOS:    "<<datEvent->fqmrdot<<endl;
    }
    cout<<"FQNSE:    "<<datEvent->fqnse<<endl;
    cout<<"FQPI0:    "<<datEvent->attribute[1]<<endl;
    if (datEvent->attribute[1]>0.){
      cout<<"  FQPI0MASS:      "<<datEvent->fqpi0mass[0]<<endl;
      cout<<"  FQPI0MOM1:      "<<datEvent->fqpi0mom1[0]<<endl;
      cout<<"  FQPI0MOM2:      "<<datEvent->fqpi0mom2[0]<<endl;
      cout<<"  FQPI0DCONV1:    "<<datEvent->fqpi0dconv1[0]<<endl;
      cout<<"  FQPI0DCONV2:    "<<datEvent->fqpi0dconv2[0]<<endl;
      cout<<"  FQPI0ANGLE:     "<<datEvent->fqpi0photangle[0]*180./TMath::Pi()<<endl;
    }
    cout<<"ENU (e):  "<<datEvent->fq1renu[0]<<endl;
    cout<<"ENU (mu): "<<datEvent->fq1renu[1]<<endl;
    if (datEvent->attribute[0]>0.){
      cout<<"E-LIKE"<<endl;
    }
    if (datEvent->attribute[0]<0.){
      cout<<"MU-LIKE"<<endl;
    }
    if (datEvent->passecut) cout<<"PASSED E CUTS"<<endl;
    else if (datEvent->passmucut) cout<<"PASSED MU CUTS"<<endl;
    else{ cout<<"FAILED"<<endl;}
    cout<<"---------------------------------------------------"<<endl;
    cout<<endl;
    cout<<endl;
    cout<<endl;
    cout<<endl;
    cout<<endl;
    cout<<endl;
    cout<<endl;
    return;
  }

  else{
   cout<<"not found"<<endl;
   return;
  }

}


////////////////////////////////////////////////////////////
void selectData::selectFV(float w1, float w2, float tw1, float tw2){

  TObjArray* list = chd->GetListOfFiles();
  TString filename;
  int Nf = list->GetEntries();
  cout<<"Nfiles: "<<Nf<<endl;
  for (int i=0; i<Nf; i++){
    filename = list->At(i)->GetTitle();
//    cout<<"Opening: "<<filename.Data()<<endl;
    TFile* f = new TFile(filename.Data());
//    cout<<"Getting Tree: "<<endl;
    TTree* tr = (TTree*)f->Get("h1");
    if (tr==NULL) continue;
    int Nev = tr->GetEntries();
//    cout<<"break"<<endl;
    fqProcessedEvent* datevent = new fqProcessedEvent(tr);
//    cout<<"break"<<endl;
    cout<<"In File: "<<filename.Data()<<" ("<<Nev<<")"<<endl;
    for (int iev=0; iev<Nev; iev++){
      tr->GetEntry(iev);
      if (datevent->passecut || datevent->passmucut){
        if (datevent->fqwall > w1 && datevent->fqwall< w2){
          if (datevent->fqtowall > tw1 && datevent->fqtowall < tw2){
            cout<<"check event: "<<iev<<"wall: "<<datevent->nev<<endl; 
          }
        }
      }
    }

    f->Close();
  }

  //
  return;
}


////////////////////////////////////////////////////////////
selectData::selectData(TChain* ch){

  chd = ch;
  datEvent = new fqProcessedEvent(chd);


}

