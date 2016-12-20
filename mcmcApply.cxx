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

  // set index of parameters
//  attIndexMom = 4;
//  attIndexRCPar = 3;
//  attIndexPID = 0;
//  attIndexPi0Par = 1;
//  attIndexPiPPar = 2;

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
// apply parameters to temporary array
void mcmcApply::applyPars(int nbin, int ncomponent, float attributeTmp[], int natt){

  float smear;
  float bias;

  for (int iatt=0; iatt<natt; iatt++){
   
    // get parameters
    smear = (float)fitPars->getAttModParameter(nbin, ncomponent, iatt, 0);
    bias = (float)fitPars->getAttModParameter(nbin, ncomponent, iatt, 1);

    // apply parameters
    attributeTmp[iatt] = smear*attributeTmp[iatt] + bias;
  }

  //
  return;
}


////////////////////////////////////////////////////////////////////////////////////////////
//
//void mcmcApply::applyPars(int nbin, int component, float &fqpidpar, float &fqmom, float &fqpi0par, float &fqpippar, float &fqrcpar);
/*
  double smear;
  double bias;


  // momentum
  // get smear parameter
  smear = fitPars->getAttModParameter(nbin, ncomponent, attIndexMom, 0);
  // get bias parameter
  bias = fitPars->getAttModParameter(nbin, ncomponent, attIndexMom, 1);
  fqmom = smear*fqmom + bias;

  // PID 
  // get smear parameter
  smear = fitPars->getAttModParameter(nbin, ncomponent, attIndexPID, 0);
  // get bias parameter
  bias = fitPars->getAttModParameter(nbin, ncomponent, attIndexPID, 1);
  fqpidpar = smear*fqpidpar + bias;


  // pip0ar 
  // get smear parameter
  smear = fitPars->getAttModParameter(nbin, ncomponent, attIndexPi0Par, 0);
  // get bias parameter
  bias = fitPars->getAttModParameter(nbin, ncomponent, attIndexPi0Par, 1);
  fqpi0par = smear*fqpi0par + bias;

  
  // pi0par
  // get smear parameter
  smear = fitPars->getAttModParameter(nbin, ncomponent, attIndexPID, 0);
  // get bias parameter
  bias = fitPars->getAttModParameter(nbin, ncomponent, attIndexPID, 1);
  fqpidpar = smear*fqpidpar + bias;

  // RCpar
  // get smear parameter
  smear = fitPars->getAttModParameter(nbin, ncomponent, attIndexPID, 0);
  // get bias parameter
  bias = fitPars->getAttModParameter(nbin, ncomponent, attIndexPID, 1);
  fqpidpar = smear*fqpidpar + bias;


  // pi0 mass
  // get smear parameter
  smear = fitPars->getAttModParameter(nbin, ncomponent, attIndexPi0Mass , 0);
  // get bias parameter
  bias = fitPars->getAttModParameter(nbin, ncomponent, attIndexPi0Mass, 1);
  fqpi0mass = smear*fqpi0mass + bias;


  // pi0 likelihood
  // get smear parameter
  smear = fitPars->getAttModParameter(nbin, ncomponent, attIndexPi0Like , 0);
  // get bias parameter
  bias = fitPars->getAttModParameter(nbin, ncomponent, attIndexPi0Like , 1);
  fqpi0like = smear*fqpi0like+ bias;

  //   
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
//    cout<<"comp: "<<event_comp<<endl;
//    cout<<"bin: "<<event_bin<<endl;
//    cout<<"bias: "<<bias<<endl;
//    cout<<"smear: "<<smear<<endl;
//    cout<<"attribute "<<iatt<<" "<<mcEvent->attribute[iatt]<<" -> ";
    mcEvent->attribute[iatt] = smear*mcEvent->attribute[iatt] + bias;
//    cout<<mcEvent->attribute[iatt]<<endl;

  }
  // otherwise, modify them all
  else{
    for (int jatt=0; jatt<fitPars->nAttributes; jatt++){
      applyPars(jatt);
    }
  }

  return;
}
*/


#endif
