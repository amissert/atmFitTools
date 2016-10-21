#ifndef MCMCAPPLY_CXX
#define MCMCAPPLY_CXX

#include "mcmcApply.h"

/////////////////////////////////////////////////////////////////////////////////////////////
// constructor
mcmcApply::mcmcApply(atmFitPars* fitpars, mcmcReader* mcmcpars, fqProcessedEvent* mcevent){

  // set pointers
  fitPars = fitpars;
  mcmcReader = mcmcreader;
  mcEvent = mcevent;

}

///////////////////////////////////////////////////////////////////////////////////////////
// set parameters to the values pointed to by mcmcpars
void mcmcApply::setFromMCMC(){

  // loop over mcmc pars and set int fitPars
  for (int ipar=0; ipar<mcmcPars->npars; ipar++){
    int atmparindex = mcmcPars->parindex[ipar];
    fitpars->setParameter(atmparindex, (double)mcmcPars->par[ipar]);
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////////////////
// apply parameters to attribute iatt of the MC event
void mcmcApply::applyPars(int iatt){

  // get this event's bin
  int event_bin = mcEvent->nbin;

  // get this event's component 
  int event_comp = mcEvent->ncomponent;

  // get smear parameter
  double smear = fitpars->getAttModParameter(event_bin, event_comp, iatt, 0);

  // get bias parameter
  double bias = fitpars->getAttModParameter(event_bin, event_comp, iatt, 1);

  // apply parameters
  mcEvent->attribute[iatt] = smear*mcEvent->attribute[iatt] + bias;

  return;
}



#endif
