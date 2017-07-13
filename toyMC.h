#ifndef TOYMC_H
#define TOYMC_H

#include "atmFitPars.h"
#include "fqProcessedEvent.h"
#include "mcmcReader.h"
#include "histoCompare.h"
#include "mcmcApply.h"
#include "modHistoArrayFV.h"
#include "histoCompare.h"
#include "eventSelectors.h"
#include "t2kSample.h"
#include "randomList.h"
#include "uniformList.h"
#include "mcLargeArray.h"
#include <exception>
#include "TH2FV.h"
#include "SKError.h"
#include "RCMap.h"


// binning for muons
const double EnuBinning[] = {0.,100.,300.,500,700,900,1100,2000,5000};
const int    EnuNBins = 8;

const double EnuBinningElectron[] = {0.,100.,300.,500,700,900,1100,2000,5000,30000};
const int    EnuNBinsElectron = 9;

// binning for muons
//const double EnuBinning[] = {0.,600.,950.,1275.,2120.,10000.};
//const double EnuBinning[] = {0.,500.,950.,1275.,2200.,20000.};
//const double EnuBinning[] = {0.,500.,1000,2000,30000.};
//const int    EnuNBins = 3;

//const double EnuBinning[] = {0.,,596.,843.,1275.,2120.};
//const int    EnuNBins = 5;

//const double EnuBinningElectron[] = {0.,400,800,1200.};
//const double EnuBinningElectron[] = {0.,1250.};
//const int    EnuNBinsElectron = 1;

//const double EnuBinning[] = {0.,100.,300.,700.,1250.,2000.,5000.,30000.};
//const int    EnuNBins = 7;

//const double EnuBinning[] = {0.,350.,800.,1250.,30000.};
//const int    EnuNBins = 4;

///////////////////////////////////////////////////
// Class to run a toy MC to apply MCMC results
class toyMC{

  public:

  toyMC();

  // vars
  TChain* chMC;
  fqProcessedEvent* mcEvent;
  TChain* chPars;
  mcmcReader* mcmcPars;
  mcmcApply* modifier;
  histoCompare* hCompare;
  atmFitPars* fitPars;
//  modHistoArrayFV* hArrFV;
  t2kSample* t2kToys;
  mcLargeArray* fastevents;
  fqcutparams cutPars;
  SKError* skErr;


  int nMCevents;
  int indexPIDPar;
  int indexPi0Par;
  int indexPiPPar;
  int indexRCPar;
//  int indexMom;

  // methods
  void setChains(TChain* chmc, TChain *chpars,int nmcevents);
  void setCompare(histoCompare* hc);
//  void setRCMap(const char* fname){rcMap = new RCMap(fname);}

//  void makeFVMapNuMu(int nmcmcpts, const char* outfile);
//  void makeFVMapNuE(int nmcmcpts,const char* outfile);
  void makeFVUncMap(int nmcmcpts, int nselection, const char* outfile, int fvbintype=1);

  void makeCombinedUncertainty(int nmcmcpts);
  void setAtmFitPars(const char* parfile);
  void fillSKErrors(int ntoys, int nbinning=1, int flgcustom=0, int effdef=0);
  void fillMarginalizedSKErr(const int ntoys, const int nmarg, int nbinning=1, int flgcustom=0, int effdef=0);
  int applyCutsToModifiedEvent(int iev,bool flgmod=true);
  int applyCoreCutsToModifiedEvent(int iev, int nclass, bool flgmod=true);
  int getEventCatagory(int iev, int nutype);
  int getRandomMCMCPoint();  

  int flgRCWeight;
};


#ifdef CINTMODE
#include "globalRandom.cxx"
#endif



#endif
