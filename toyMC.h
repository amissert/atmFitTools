#ifndef TOYMC_H
#define TOYMC_H

#include "atmFitPars.h"
#include "fqProcessedEvent.h"
#include "mcmcReader.h"
#include "histoCompare.h"
#include "mcmcApply.h"
#include "modHistoArray.h"
#include "histoCompare.h"
#include "eventSelector.h"
#include "TH2FV.h"


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
  modHistoArray* hArr;

  // methods
  void setChains(TChain* chmc, TChain *chpars);
  void setCompare(histoCompare* hc);
  void fillArrayDirect(int isamp, int ibin, int iatt, int npts);
  void testToy(int nmcmcpts);
  int getRandomMCMCPoint();  

};


#ifdef CINTMODE
#include "globalRandom.cxx"
#endif



#endif
