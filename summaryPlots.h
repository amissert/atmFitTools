#ifndef SUMMARYPLOTS_H
#define SUMMARYPLOTS_H

#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "mcLargeArray.h"


#define NATTS 5
#define NTOYPOINTS 50
#define NCATS 5

//////////////////////////////////////////////////////////////////
// Container for useful plots to summarize event selection results 
class summaryPlots{
  public:

  // constructor
  summaryPlots(const char* name = "");

  // vars
  TString nameTag; 


//  TH1D* pltToySpectrum[NTOYPOINTS];
//  TH1D* pltToyPower[NTOYPOINTS];
//  TH1D* pltToySyst[NTOYPOINTS];

  // nuMu hisots
  TH1D* pltEnuMu; //< rc spectrum
  TH1D* pltEnuMuCat[NCATS]; //< rc spectrum
  TH1D* pltEnuMuSg;
  TH1D* pltEnuMuBg;
  TH1D* pltPowerMu; //< oscillation power
  TH1D* pltPassMu; //< passed nominal cuts?

  // nuE histos
  TH1D* pltEnuE; //< rc spectrum
  TH1D* pltEnuESg;
  TH1D* pltEnuEBg;
  TH1D* pltEnuECat[NCATS]; //< rc spectrum
  TH1D* pltPowerE; //< oscillation power
  TH1D* pltPassE; //< passed nominal cuts?

  // nuE 1RPi histos
  TH1D* pltEnuE1RPi; //< rc spectrum
  TH1D* pltEnuE1RPiSg;
  TH1D* pltEnuE1RPiBg;
  TH1D* pltEnuE1RPiCat[NCATS]; //< rc spectrum
  TH1D* pltPowerE1RPi; //< oscillation power
  TH1D* pltPassE1RPi; //< passed nominal cuts?

  // other histos
  TH1D* pltSyst; //< additional systematic uncertainty
  TH1D* pltAtt[NATTS]; //< fq cut attributes
  TH2D* pltPower2D; //< osc power vs. true energy;
  mcLargeArray* fastevents;

  // initialize all histos
  void Init();

  // initialize toy array
//  void InitToys(TH1D* hseed);

  // set pointer to large mc array
  void setLargeArray(mcLargeArray* lgarr);

  // fill all histos from large array 
  void fillAllFromArray(int iev, float pow=-1, float sys=-1);
  
  // get the event catagory
  int GetCatagory(int iev, int wantnutype);

  // clear all the histos
  void clearHistos();

  // set poisson errors
  void setPoissonErrors(TH1D* hh);

  // find the depth of the oscillation dip for numu
  float calcDipDepth();

  // get the signal and background histograms
//  void calcSignalAndBg(
};





#ifdef CINTMODE
#include "summaryPlots.cxx"
#endif


#endif












