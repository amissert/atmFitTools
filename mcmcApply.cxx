#ifndef MCMCAPPLY_CXX
#define MCMCAPPLY_CXX

#include "mcmcApply.h"

/////////////////////////////////////////////////////////////////////////////////////////////
// constructor
mcmcApply::mcmcApply(atmFitPars* fitpars, mcmcReader* mcmcpars){

  // set pointers
  fitPars = fitpars;
  mcmcPars = mcmcpars;

  // attribute indicies
  indexPIDPar = 0;
  indexPi0Par = 1;
  indexPiPPar = 2;
  indexRCPar  = 3;

}

///////////////////////////////////////////////////////////////////////////////////////////
// set parameters to the values pointed to by mcmcpars
void mcmcApply::setFromMCMC(){

  // loop over mcmc pars and set int fitPars
  for (int ipar=0; ipar<mcmcPars->npars; ipar++){
   
    // get the index of the ipar-th parameter in the atmfitpars array
    int atmparindex = mcmcPars->parindex[ipar];

    // talk about it
    cout<<"set par "<<atmparindex<<" "<<" <- "<<ipar<<" "<<mcmcPars->par[ipar]<<endl;

    // change the atmfitpars array
    fitPars->setParameter(atmparindex, (double)mcmcPars->par[ipar]);
  }

  return;
}


////////////////////////////////////////////////////////////////////////////////////////////
// apply parameters to temporary array
void mcmcApply::applyPars(int nbin, int ncomponent, float attributeTmp[], int natt){

  float smear;
  float bias;
  
//  cout<<"----"<<endl;
  for (int iatt=0; iatt<natt; iatt++){
   
    // get parameters
    smear = (float)fitPars->getAttModParameter(nbin, ncomponent, iatt, 0);
    bias = (float)fitPars->getAttModParameter(nbin, ncomponent, iatt, 1);
//    cout<<"bias: "<<bias<<endl;
//    cout<<"smear: "<<smear<<endl;
    // apply parameters
//    cout<<"att "<<iatt<<": "<<attributeTmp[iatt]<<" -> ";
    attributeTmp[iatt] = smear*attributeTmp[iatt] + bias;
    cout<<attributeTmp[iatt]<<endl;
  }

  //
  return;
}


/////////////////////////////////////////////////////////////////
// apply the cuts to a modified event and see if it passes
int mcmcApply::applyCutsToModifiedEvent(int iev, mcLargeArray* fastevents, bool modflg){

  // fill tmp array with "nominal" MC values
  const int natt = 4;
  float attributesTmp[natt];
  for (int iatt=0; iatt<natt; iatt++){
    attributesTmp[iatt] = fastevents->vattribute[iev][iatt];   
  }
 
  // modify tmp array by applying the histogram shape parameters
  if (modflg){ applyPars(fastevents->vbin[iev],
               fastevents->vcomponent[iev],
               attributesTmp,
               natt);}

  // structure for T2K cuts
  fqcutparams cutPars;

  // fill cut parameter structure using modified attributes
  cutPars.fqpid = attributesTmp[indexPIDPar];
  cutPars.fqpi0par = attributesTmp[indexPi0Par];
  cutPars.fqpippar = attributesTmp[indexPiPPar];
//  cutPars.fqrcpar = attributesTmp[indexRCPar];

  // fill cut parameter structure using modified attributes
//  cutPars.fqpid = attributesTmp[0];
//  cutPars.fqpi0par = attributesTmp[1];
//  cutPars.fqpippar = attributesTmp[2];
//  cutPars.fqrcpar = attributesTmp[3];

  // other cut pars that are not modified
  cutPars.fqmome = fastevents->vfqmumom[iev];
  cutPars.fqmommu = fastevents->vfqemom[iev];
  cutPars.nhitac = fastevents->vnhitac[iev];
  cutPars.fqnsubev = fastevents->vfqnsubev[iev];
  cutPars.fqenue = fastevents->vfqenue[iev];
  cutPars.fqenumu = fastevents->vfqenumu[iev];
  cutPars.fqrcpar = fastevents->vfqrcpar[iev];
  cutPars.fqnring = fastevents->vfqnring[iev];

  // see if it passes cuts
  int passnue = selectNuE(cutPars);
  int passnumu = selectNuMu(cutPars);
  
  //
  if (passnue>0) return 1;
  if (passnumu>0) return 2;

  //
  return 0;
  
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
