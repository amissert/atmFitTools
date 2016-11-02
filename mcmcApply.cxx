#ifndef MCMCAPPLY_CXX
#define MCMCAPPLY_CXX

#include "mcmcApply.h"

/////////////////////////////////////////////////////////////////////////////////////////////
// constructor
mcmcApply::mcmcApply(atmFitPars* fitpars, mcmcReader* mcmcpars, fqProcessedEvent* mcevent){

  // set pointers
  fitPars = fitpars;
  mcmcPars = mcmcpars;
  mcEvent = mcevent;

}

///////////////////////////////////////////////////////////////////////////////////////////
// set parameters to the values pointed to by mcmcpars
void mcmcApply::setFromMCMC(){

  // loop over mcmc pars and set int fitPars
  for (int ipar=0; ipar<mcmcPars->npars; ipar++){
    int atmparindex = mcmcPars->parindex[ipar];
    cout<<"set par "<<atmparindex<<" "<<" <- "<<ipar<<" "<<mcmcPars->par[ipar]<<endl;
    fitPars->setParameter(atmparindex, (double)mcmcPars->par[ipar]);
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

  // modify single parameter if specified
  if (iatt>=0){

    // get smear parameter
    double smear = fitPars->getAttModParameter(event_bin, event_comp, iatt, 0);

    // get bias parameter
    double bias = fitPars->getAttModParameter(event_bin, event_comp, iatt, 1);

    // apply parameters
    cout<<"comp: "<<event_comp<<endl;
    cout<<"bin: "<<event_bin<<endl;
    cout<<"bias: "<<bias<<endl;
    cout<<"smear: "<<smear<<endl;
    cout<<"attribute "<<iatt<<" "<<mcEvent->attribute[iatt]<<" -> ";
    mcEvent->attribute[iatt] = smear*mcEvent->attribute[iatt] + bias;
    cout<<mcEvent->attribute[iatt]<<endl;

  }
  // otherwise, modify them all
  else{
    for (int jatt=0; jatt<fitPars->nAttributes; jatt++){
      applyPars(jatt);
    }
  }

  return;
}



#endif
