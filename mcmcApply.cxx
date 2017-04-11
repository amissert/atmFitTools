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

  // flags
  flgApplyXSecPar=false;
  flgApplyFluxPar=false;
  flgApplyNormPar=true;
  flgApplyTest=false;
  for (int iatt=0; iatt<NATTMAX; iatt++){
    flgApplyAttSmearPar[iatt] = true;
    flgApplyAttBiasPar[iatt]  = true;
  }
  flgGlobalUseBestPars=false;
  for (int ipar=0; ipar<NMCMCPARS; ipar++){
    flgUseBestPar[ipar]=false;
  }

  findBestFitPars();

}



///////////////////////////////////////////////////////////////////////////////////////////
// find the best fit pars (used to calculate SK errrors)
void mcmcApply::findBestFitPars(){

   // setup array
   const int NFitPars = fitPars->nTotPars;
   const int NMCMCSamples = mcmcPars->fChain->GetEntries();
   float par_array[NFitPars][NMCMCSamples];

   // fill array
   cout<<"mcmcApply::findBestFitPars(): finding the best fit parameters..."<<endl;
   for (int isamp=0; isamp<NMCMCSamples; isamp++){
     mcmcPars->GetEntry(isamp);
     for (int jpar =0; jpar<NFitPars; jpar++){
       par_array[jpar][isamp] = mcmcPars->par[jpar]; 
     }
   }

   // calculate means
   for (int ipar =0; ipar<NFitPars; ipar++){
      bestFitPar[ipar]  = arraymean(par_array[ipar],NMCMCSamples);
   }

   return;
}



///////////////////////////////////////////////////////////////////////////////////////////
void mcmcApply::setUseBestFitSystPars(bool value){

  // loop ofer sys pars and set flag to value
  for (int i=0; i<(fitPars->nSysPars - fitPars->nNormPars); i++){
    int index = fitPars->getSysParIndex(i); 
    cout<<"mcmcApply::setUseBestFitSystPars: Fixing par"<<index<<endl;
    flgUseBestPar[index] = value; 
  }
  
  //
  return;
}



///////////////////////////////////////////////////////////////////////////////////////////
void mcmcApply::setUseBestFitNormPars(bool value){

  // loop over norm pars and set flag to value
  for ( int isamp=0; isamp<fitPars->nSamples; isamp++ ) {
    for ( int jbin=0; jbin<fitPars->nBins; jbin++ ) {
      int index = fitPars->getNormParIndex(isamp,jbin);
      cout<<"mcmcApply::setUseBestFitNormPars: Fixing par"<<index<<endl;
      flgUseBestPar[index]=value;
    }
  }

  //
  return;
}




///////////////////////////////////////////////////////////////////////////////////////////
void mcmcApply::setBestPars(){

  // loop over mcmc pars and set int fitPars
  for (int ipar=0; ipar<mcmcPars->npars; ipar++){

    // replace with best fit if using best fit for some parameters
    if (flgGlobalUseBestPars){
      if (flgUseBestPar[ipar]){
        fitPars->setParameter(ipar, bestFitPar[ipar]);
      }
    }
  }

  return;
}



///////////////////////////////////////////////////////////////////////////////////////////
// set parameters to the values pointed to by mcmcpars
void mcmcApply::setAlphaFromMCMC(){

  int n_beta = fitPars->nTotPars - fitPars->nSysPars;

  // loop over mcmc pars and set int fitPars
  for (int ipar=0; ipar<mcmcPars->npars; ipar++){
  
    // get the index of the ipar-th parameter in the atmfitpars array
    int atmparindex = mcmcPars->parindex[ipar];

    // set if it's alpha
    if (atmparindex >= n_beta){

      // talk about it
//      cout<<"set par "<<atmparindex<<" "<<" <- "<<ipar<<" "<<mcmcPars->par[ipar]<<endl;

      // change the atmfitpars array
      fitPars->setParameter(atmparindex, (double)mcmcPars->par[ipar]);

    }

  }

  return;
}



///////////////////////////////////////////////////////////////////////////////////////////
// set parameters to the values pointed to by mcmcpars
void mcmcApply::setBetaFromMCMC(){

  //
  int n_beta = fitPars->nTotPars - fitPars->nSysPars;

  // loop over mcmc pars and set int fitPars
  for (int ipar=0; ipar<mcmcPars->npars; ipar++){
  
    // get the index of the ipar-th parameter in the atmfitpars array
    int atmparindex = mcmcPars->parindex[ipar];

    // set if it's alpha
    if (atmparindex < n_beta){

      // talk about it
//      cout<<"set par "<<atmparindex<<" "<<" <- "<<ipar<<" "<<mcmcPars->par[ipar]<<endl;

      // change the atmfitpars array
      fitPars->setParameter(atmparindex, (double)mcmcPars->par[ipar]);

    }

  }


  return;
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
// apply parameters in temporary array
void mcmcApply::applyPars(int nbin, int ncomponent, float attributeTmp[], int natt){

  
  for (int iatt=0; iatt<natt; iatt++){

    float smear=1.0;
    float bias=0.0;

    // get parameters
    if (flgApplyAttSmearPar[iatt]){
      smear = (float)fitPars->getAttModParameter(nbin, ncomponent, iatt, 0);
    }
    if (flgApplyAttBiasPar[iatt]){
      bias = (float)fitPars->getAttModParameter(nbin, ncomponent, iatt, 1);
    }

    attributeTmp[iatt] = smear*attributeTmp[iatt] + bias;


  }

  //
  return;
}



////////////////////////////////////////////////////////////////////////////////////////////
// get total weight from atm flux pars
float mcmcApply::getFluxWeight(float enu, int nutype){

  float ww = 1.0;
  
  if (enu<1000.){
    ww*= fitPars->getSysParameter(14);
  }

  else if (enu>=1000.){
    ww*= fitPars->getSysParameter(15);
  }

  if (nutype==14){
    ww*=fitPars->getSysParameter(18);
  }
 
//  cout<<"enu: "<<enu;
//  cout<<" nutype: "<<nutype;
//  cout<<" flxweight: "<<ww<<endl;;

  return ww;
}



////////////////////////////////////////////////////////////////////////////////////////////
// get total weight from xsec pars
float mcmcApply::getXsecWeight(int mode, float Enu){

  // weights from tn186!! ////////////////////////////
  float ww = 1.0; 
  if (mode<10){
      if (Enu<190.) ww*=fitPars->getSysParameter(0);
      if ((Enu>190.)&&(Enu<240.)) ww*=fitPars->getSysParameter(1);
      if ((Enu>240.)&&(Enu<294.)) ww*=fitPars->getSysParameter(2);
      if ((Enu>294.)&&(Enu<333.)) ww*=fitPars->getSysParameter(3);
      if ((Enu>333.)&&(Enu<390.)) ww*=fitPars->getSysParameter(4);
      if ((Enu>390.)&&(Enu<440.)) ww*=fitPars->getSysParameter(5);
      if ((Enu>440.)&&(Enu<487.)) ww*=fitPars->getSysParameter(6);
      if ((Enu>487.)&&(Enu<590.)) ww*=fitPars->getSysParameter(7);
      if ((Enu>590.)&&(Enu<690.)) ww*=fitPars->getSysParameter(8);
      if ((Enu>690.)&&(Enu<786.)) ww*=fitPars->getSysParameter(9);
      if ((Enu>786.)&&(Enu<896.)) ww*=fitPars->getSysParameter(10);
      if ((Enu>896.)&&(Enu<994.)) ww*=fitPars->getSysParameter(11);
      if ((Enu>994.)&&(Enu<2000.)) ww*=fitPars->getSysParameter(12);
      if (Enu>2000.) ww*=fitPars->getSysParameter(13);
  }
  else if (mode<30){
    ww *= fitPars->getSysParameter(16);
  }
  else if (mode>=30){
    ww *= fitPars->getSysParameter(17);
  }

  //
  return ww;

}



////////////////////////////////////////////////////////////////////////////////////////////
// can apply a test custom weight here
float mcmcApply::getTestWeight(int nbin, int nsamp, int nmode, float enutrue, int nutype){

  // start fresh
  float ww = 1.0;

  if (nutype==14){
    ww*= fitPars->getSysParameter(14);
  }

  return ww;
}



////////////////////////////////////////////////////////////////////////////////////////////
// apply parameters to temporary array
float mcmcApply::getEvtWeight(int nbin, int nsamp, int nmode, float enutrue, int nutype){

  // start fresh
  float ww = 1.0;

  // apply normalization
  if (flgApplyNormPar) ww *= (float)fitPars->getNormParameter(nsamp,nbin);

  // apply xsec
  if (flgApplyXSecPar) ww *= getXsecWeight(nmode, enutrue);

  // apply flux
  if (flgApplyFluxPar) ww *= getFluxWeight(enutrue,nutype);

  if (flgApplyTest) ww *= getTestWeight(nbin,nsamp,nmode,enutrue,nutype);
  return ww;
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
  cutPars.fqrcpar = attributesTmp[indexRCPar];

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



#endif
