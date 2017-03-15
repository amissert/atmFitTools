
#ifndef SKERROR_H
#define SKERROR_H



#include <iostream>

#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TChain.h"
#include "TMath.h"
#include "TString.h"
#include "stats.C"
#include "TLine.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLatex.h"


const int NCLASSES = 50;
const int NTOYS    = 2000;
const int NLINES= 3;

using namespace std;


///////////////////////////////////////////////////////////////////
//
class SKError{

  public:

  // constructor
  SKError(int ntoy = 100);
  int Nclass;
  int Ntoys;
 
  // arrays for numbers of events
  double Nevents[NCLASSES][NTOYS];
  double NeventsTotal[NCLASSES][NTOYS];
//  double Efficiency[NCLASSES][NTOYS];
  double DelEfficiency[NCLASSES][NTOYS];
  double DelEffShiftError[NCLASSES];
  double absDelEffShiftError[NCLASSES];
  double DelEffFitError[NCLASSES];


  // evis binning
  TH1D* hEvisNuECCQE;
  TH1D* hEvisNuECCOth;
  TH1D* hEvisNuMuCCQE;
  TH1D* hEvisNuMuCCOth;

  // totals
  TH1D* hEvisNuECCQETot;
  TH1D* hEvisNuECCOthTot;
  TH1D* hEvisNuMuCCQETot;
  TH1D* hEvisNuMuCCOthTot;

  // for drawing lines
  TLine* lineHorz[NLINES];
  TLine* lineVert[NLINES];
  TLatex* labelVert[NCLASSES];
  TLatex* labelHorz[NCLASSES];
  TLatex* sectorLabelHorz[4];
  TLatex* sectorLabelVert[4];
  TLatex* nuLabelVert[4];
  TLatex* nuLabelHorz[4];
  int    lineVal[NLINES];  

  TH2D* hCor;
  TH2D* hCov;
  TVectorD* vShiftErrors;
  TH1D* hDiagonalErrors;

  // histogram of all numbers of events
  TH1D* hSlice;
  TGraph* gScat;
  TH1D* hdist;
  TLine* distMean;
  TLine* zeroValue;


  // initialize histograms
  void initHistos(int ibinning=0);
  void zeroArrays();
  void resetHistos();

  // draw a particular toy
  void drawSlice(int ntoy);
  void drawSliceTot(int ntoy);
  void drawSliceEff(int ntoy);
  void drawAll();
  void drawAllEff();
  void drawDist(int nclass);
  void drawEffDist(int nclass);
  void drawCor();
  void drawCov();


  // for classifying and filling
  int getClassMC(int nutype, int mode, int component, double evis, int nsubev, double towall, double wall);

  // fill an event in the histograms
  void addEvent(int nclass, double evis, double weight, bool flgtotal); 

  // save histo contents into arrays
  void addToy(int ntoy);

  // calculate correlation and covariance
  void calcCov();

  // calculate correlation and covariance using epsilon
  void calcCovEff();

  // calculate correlation and covariance using epsilon
  void calcCovDelEff();

  // calculate efficiency based on index of total event numbers
  double calcEff(int nclass, int ntoy);

  // calculate efficiency based on index of total event numbers
  double calcDelEff(int nclass, int ntoy);

  // calculate all effeciencies
  void calcAllEff(int ntoy);;

  // calculate all effeciencies
  void calcAllDelEff(int ntoy);;

  // draw scatterplot
  void drawScatter(int iclass, int jclass);

  //
  double calcShiftError(int iclass){
    cout<<"Shift error for class: "<<iclass<<endl;
    return arraymeanD(DelEfficiency[iclass],Ntoys);
  }

  //
  double calcFitError(int iclass){
    cout<<"Fit error for class: "<<iclass<<endl;
    return TMath::Sqrt(arrayvarD(DelEfficiency[iclass],Ntoys));
  }

  void calcErrors(){
    for (int iclass=0; iclass<Nclass; iclass++){
      DelEffShiftError[iclass] = calcShiftError(iclass);
      absDelEffShiftError[iclass] = TMath::Abs(calcShiftError(iclass));
      DelEffFitError[iclass] = calcFitError(iclass);
    }
    vShiftErrors = new TVectorD(Nclass,DelEffShiftError);
    return;
  }

  // 
  void printErrors(){
    for (int iclass=0; iclass<Nclass; iclass++){
      cout<<"--- Class "<<iclass<<"  ---"<<endl;
      cout<<"  Fit: "<<DelEffFitError[iclass]*100.<<"%";
      cout<<"  Shift: "<<DelEffShiftError[iclass]*100.<<"%";
      double toterr = TMath::Sqrt(DelEffFitError[iclass]*DelEffFitError[iclass]
                                 +DelEffShiftError[iclass]*DelEffShiftError[iclass]);
      cout<<"  Total: "<<toterr*100.<<"%"<<endl;
    }
  }

  void saveErrors(const char* filename){
    TFile* fout = new TFile(filename,"RECREATE");
    hCov->Write();
    hCor->Write();
    vShiftErrors->Write("vshift");
    cout<<"writing: "<<filename<<endl;
    fout->Close();
    return;
  }

  void printEffDist(const char* plotdir){
     TCanvas* cc = new TCanvas("cc","cc",700,600);
     cc->GetPad(0)->SetLeftMargin(0.12);
     cc->GetPad(0)->SetRightMargin(0.12);
     cc->GetPad(0)->SetBottomMargin(0.15);
     TString pdir = plotdir;
     for (int i=0; i<Nclass; i++){
       TString pname = pdir.Data();
       pname.Append(Form("throw_dist_class_%d.png",i));
       drawEffDist(i);
       cc->Print(pname.Data());
     }
     return;
  }

};



#ifdef CINTMODE
#include "SKError.cxx"
#endif










#endif
