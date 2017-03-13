#ifndef TOYMC_CXX
#define TOYMC_CXX

#include "toyMC.h"

using namespace std;


/////////////////////////////////////////////////////////////////
// apply the cuts to a modified event and see if it passes
//   returns 1 -> passed electron selection
//   returns 2 -> passed muon selection
int toyMC::applyCutsToModifiedEvent(int iev, bool flgmod){

  // fill tmp array with "nominal" MC values
  const int natt = 4;
  float attributesTmp[natt];
  for (int iatt=0; iatt<natt; iatt++){
    attributesTmp[iatt] = fastevents->vattribute[iev][iatt];   
  }

  if (indexPIDPar>=0) cutPars.fqpid = attributesTmp[indexPIDPar];
  if (indexPi0Par>=0) cutPars.fqpi0par = attributesTmp[indexPi0Par];
  if (indexPiPPar>=0) cutPars.fqpippar = attributesTmp[indexPiPPar];
  if (indexRCPar>=0) cutPars.fqrcpar = attributesTmp[indexRCPar];

  // modify tmp array by applying the histogram shape parameters
  if (flgmod){
    // modify these cut values by the shape parameters
    modifier->applyPars(fastevents->vbin[iev],
                        fastevents->vcomponent[iev],
                        attributesTmp,
                        natt);

   // fill cut parameter structure using modified attributes
   if (indexPIDPar>=0) cutPars.fqpid = attributesTmp[indexPIDPar];
   if (indexPi0Par>=0) cutPars.fqpi0par = attributesTmp[indexPi0Par];
   if (indexPiPPar>=0) cutPars.fqpippar = attributesTmp[indexPiPPar];
   if (indexRCPar>=0) cutPars.fqrcpar = attributesTmp[indexRCPar];

  }

  // other cut pars that are not modified
  cutPars.fqmome = fastevents->vfqmumom[iev];
  cutPars.fqmommu = fastevents->vfqemom[iev];
  cutPars.nhitac = fastevents->vnhitac[iev];
  cutPars.fqnsubev = fastevents->vfqnsubev[iev];
  cutPars.fqenue = fastevents->vfqenue[iev];
  cutPars.fqenumu = fastevents->vfqenumu[iev];
  cutPars.fqnring = fastevents->vfqnring[iev];

  // see if it passes cuts
  int passnue = selectNuE(cutPars);
  int passnumu = selectNuMu(cutPars);
  int passnue1rpi = selectNuE1Rpi(cutPars);

  //
  if (passnue>0) return 1;
  if (passnumu>0) return 2;
  if (passnue1rpi>0) return 3;

  //
  return 0;
  
}



/////////////////////////////////////////////////////////////////
// Same as above but applies "core" cuts to modified event (a la
// TN 186)
//   returns:
//     -1 -> Neither core nor tail
//      0 -> Tail
//      1 -> Core
/////////////////////////////////////////////////////////////////
int toyMC::applyCoreCutsToModifiedEvent(int iev, int nclass, bool flgmod){

  // fill tmp array with "nominal" MC values
  const int natt = 4;
  float attributesTmp[natt];
  for (int iatt=0; iatt<natt; iatt++){
    attributesTmp[iatt] = fastevents->vattribute[iev][iatt];   
  }

  // also set cut parameters structure
  if (indexPIDPar>=0) cutPars.fqpid = attributesTmp[indexPIDPar];
  if (indexPi0Par>=0) cutPars.fqpi0par = attributesTmp[indexPi0Par];
  if (indexPiPPar>=0) cutPars.fqpippar = attributesTmp[indexPiPPar];
  if (indexRCPar>=0) cutPars.fqrcpar = attributesTmp[indexRCPar];

  // modify tmp array by applying the histogram shape parameters
  if (flgmod){
    modifier->applyPars(fastevents->vbin[iev],
                        fastevents->vcomponent[iev],
                        attributesTmp,
                        natt);
  }

  // re-fill cut parameter structure using modified attributes
  if (indexPIDPar>=0) cutPars.fqpid = attributesTmp[indexPIDPar];
  if (indexPi0Par>=0) cutPars.fqpi0par = attributesTmp[indexPi0Par];
  if (indexPiPPar>=0) cutPars.fqpippar = attributesTmp[indexPiPPar];
  if (indexRCPar>=0) cutPars.fqrcpar = attributesTmp[indexRCPar];
  cutPars.fqnring = fastevents->vfqnring[iev];

  // other cut pars that are not modified
  cutPars.fqmome = fastevents->vfqmumom[iev];
  cutPars.fqmommu = fastevents->vfqemom[iev];
  cutPars.nhitac = fastevents->vnhitac[iev];
  cutPars.fqnsubev = fastevents->vfqnsubev[iev];
  cutPars.fqenue = fastevents->vfqenue[iev];
  cutPars.fqenumu = fastevents->vfqenumu[iev];
  cutPars.fqnring = fastevents->vfqnring[iev];

  // see if it passes core cuts
  // classes: 1 -> nu e CCQE
  //          2 -> nu mu CCQE
  //          3 -> nu e CCOth
  //          4 -> nu mu CCOth
  //          5 -> NC pi0

  if (nclass==1){ //< require single ring electron and CCQE
    if (cutPars.fqpid>=0. && cutPars.fqpi0par<=0. && cutPars.fqrcpar<=0.){//< e-like, not pi0, 1R-like
      return 1;
    }
    else{
      return 0;
    }
  }


  // muon CCQE
  if (nclass==2){
    if (cutPars.fqpid<=0 && cutPars.fqpippar<=0. && cutPars.fqrcpar<=0.){
      return 1; // mu-lik and not pip and 1R-like so is core
    }
    else{
      return 0;
    }
  }


  // electron CCOth
  if (nclass==3){
    if (cutPars.fqpid>=0 && cutPars.fqpi0par<=0. && cutPars.fqrcpar<=0.){
      return 1; // e-like and not pi0 and 1R-like
    }
    else{
      return 0;
    }
  }


  // muon CCOth
  if (nclass==4){
    if (cutPars.fqpid<=0 && cutPars.fqpippar<=0. && cutPars.fqrcpar<=0.){
      return 1; //< mu-lik and not pip and 1R-like
    }
    else{
      return 0.;
    }
  }


  // NC
  if (nclass==5){
    if (cutPars.fqpid>=0 && cutPars.fqpi0par<=0. && cutPars.fqrcpar<=0.){
      return 1; //< e-like and not pi0-like
    }
    else{
      return 0.;
    }
  }

  //
  return -1;
   
  
}



/////////////////////////////////////////////////////////////////
toyMC::toyMC(){


}

/////////////////////////////////////////////////////////////////
// set the pointers to the mc events and the post-mcmc parameters
void toyMC::setChains(TChain* chmc, TChain* chpars, int nmcevents){

  chMC = chmc;
  chPars = chpars;

  mcEvent = new fqProcessedEvent(chMC);
  fastevents = new mcLargeArray(chMC,nmcevents);
  nMCevents = nmcevents;

  mcmcPars = new mcmcReader(chPars);

  return;
}

////////////////////////////////////////////////////////////////
//  read parameters from a random point in the MCMC cloud
int toyMC::getRandomMCMCPoint(){
  int nmax = chPars->GetEntries();
  int randpoint = randy->Integer(nmax);
  cout<<"getting mcmc point: "<<randpoint<<endl;
  chPars->GetEntry(randpoint);
  return randpoint;
}



/*
////////////////////////////////////////////////////////////////
// make map of uncertainties for nu mu events
void toyMC::makeFVMapNuMu(int nmcmcpts, int mcevents){

  // make array of histos
  cout<<"Initializing array of histograms..."<<endl;
  TH1D* hE = new TH1D("hE","hE",20,0,2000);
  TH2FV* hfv = new TH2FV("hfv",1);
  // array of nu energy histograms
  hArrFV = new modHistoArrayFV(hE,hfv,nmcmcpts);

  // get list of mc events
  int nevmax = chMC->GetEntries();
  if (mcevents<nevmax) nevmax = mcevents;

  // get list of points in mcmc parameter space
  cout<<"Making list of MCMC points"<<endl;
  randomList* mcmclist = new randomList(nmcmcpts,chPars->GetEntries(),nmcmcpts);

  // fill array of T2K MC
//  mcLargeArray* fastevents = new mcLargeArray(chmc,nevmax);

  // loop over mcmc points
  for (int i=0; i<nmcmcpts; i++){

    // read in parameters
    chPars->GetEntry(mcmclist->getAt(i));

    // modify attributes using thes parameters
    modifier->setFromMCMC();

    // loop over T2K MC events
    for (int iev=0; iev<nevmax; iev++){

      // read in event
      if ((iev%5000)==0) cout<<"toyMC: getting event: "<<iev<<endl;
      chMC->GetEntry(iev);
      
      // apply the mcmc parameters
      modifier->applyPars();

      // see if event passes cuts
      float lmom = (float)mcEvent->attribute[3]; // best momentum
      float emom = (float)mcEvent->fq1rmom[0][1];
      float lpid = (float)mcEvent->attribute[0]; // PID likelihood ratio
      float enu  = (float)calcENu(2,lmom,mcEvent->fq1rdir[0][2][0],mcEvent->fq1rdir[0][2][1],mcEvent->fq1rdir[0][2][2]);
      int ipass = selectNuMu(mcEvent->nhitac,
                             mcEvent->fqnse,
                             enu,
                             lmom,
                             emom,
                             lpid,
                             mcEvent->fqmrnring[0]);

      // if it passes fill histos
      if (ipass!=1) continue;
      int fvbin = hArrFV->hFV[i]->Fill(mcEvent->fq1rtowall[0][2],mcEvent->fq1rwall[0][2],mcEvent->evtweight) - 1;
      if (fvbin>=0) hArrFV->getHistogram(i,fvbin)->Fill(enu, mcEvent->evtweight);
    }
  }
  return;

}
*/

////////////////////////////////////////////////////////////////
// get thec combined uncertainties for all events
void toyMC::makeCombinedUncertainty(int nmcmcpts){

  // setup containter for t2k sample
  t2kToys = new t2kSample("_toymc",1,1);

  // get list of mc events
  int nevmax = chMC->GetEntries();
  if (nMCevents>nevmax) nMCevents = nevmax;

  // get list of points in mcmc parameter space
  cout<<"Making list of MCMC points"<<endl;
  randomList* mcmclist = new randomList(nmcmcpts,chPars->GetEntries(),nmcmcpts);

  // loop over mcmc points
  for (int i=0; i<nmcmcpts; i++){

    // read in shape parameters
    cout<<"getting event"<<mcmclist->getAt(i)<<endl;
    chPars->GetEntry(mcmclist->getAt(i));

    // modify attributes using thes parameters
    modifier->setFromMCMC();

    // loop over T2K MC events
    for (int iev=0; iev<nMCevents; iev++){
     
      // apply cuts
      int ipass = applyCutsToModifiedEvent(iev);
      if (ipass==0) continue;

      // fill histos
      if (ipass==1) t2kToys->hEnuElectron->Fill(fastevents->vfqenue[iev], fastevents->vweight[iev]);
      if (ipass==2) t2kToys->hEnuMuon->Fill(fastevents->vfqenumu[iev], fastevents->vweight[iev]);
    }
    t2kToys->finishToyRun();
  }
  
  t2kToys->calcUncertainties();

  return;

}



////////////////////////////////////////////////////////////////
// fill the SKError class
void toyMC::fillSKErrors(int ntoys,int nbinning, int flgcustom){

  // make error container
  skErr = new SKError(ntoys);
  skErr->initHistos(nbinning);


  if (!flgcustom){
    modifier->flgApplyXSecPar = true;
    modifier->flgApplyFluxPar = true;
    modifier->flgApplyNormPar= true;
  }

  // get list of random MCMC points from the MCMC point file
//  cout<<"Making list of MCMC points"<<endl;
//  randomList* mcmclist = new randomList(ntoys,chPars->GetEntries(),ntoys);
  // sample uniforly instead
  int nmcmcmax = chPars->GetEntries();
  int nskip = nmcmcmax/ntoys;

  // determine max events
  int nevmax = chMC->GetEntries();
  if (nMCevents>nevmax) nMCevents = nevmax;

  // loop over mcmc points
  int i=0;
  while (true){

    // read in fit parameters
    int ientry = nskip*i;
    cout<<"getting event"<<ientry<<endl;
    if (ientry>=nmcmcmax) break; 
    chPars->GetEntry(ientry);

    // modify attributes using thes parameters
    if (i>0) modifier->setFromMCMC();

    // reset histogram contents
    skErr->resetHistos();

    // loop over T2K MC events
    for (int iev=0; iev<nMCevents; iev++){

      // get true MC event class for event "iev"
      int nclass = skErr->getClassMC(fastevents->vnutype[iev],
                                     fastevents->vmode[iev],
                                     fastevents->vcomponent[iev],
                                     fastevents->vfqemom[iev],
                                     fastevents->vfqnsubev[iev],
                                     fastevents->vfqtowall[iev],
                                     fastevents->vfqwall[iev]);
    
      // passes core selection?
      int iscore = applyCoreCutsToModifiedEvent(iev,nclass,true);

      // if its core or tail fill histos
      if (iscore<0) continue; //< skip unclassified events

      // get the new event weight
      float   ww = fastevents->vweight[iev]
                   *modifier->getEvtWeight( fastevents->vbin[iev], fastevents->vsample[iev]
                                           ,fastevents->vmode[iev], fastevents->vpmomv[iev]
                                           ,fastevents->vnutype[iev]);
     

      // fill total histos
      // the "true" flag means we add it to the "total" evis histogram for this class
      skErr->addEvent(nclass,fastevents->vfqemom[iev],ww,true);

      if (iscore==1){
        // fill core histos
      // the "true" flag means we add it to the "core" evis histogram for this class
        skErr->addEvent(nclass,fastevents->vfqemom[iev],ww,false);
      }
     
    }

    // save toy histo contents
    skErr->addToy(i);
    i++;;
  }

  skErr->calcCovEff();
  skErr->calcErrors();
  cout<<"filled "<<i<<" toys!"<<endl;
  //
  return;
}


/*
////////////////////////////////////////////////////////////////
// make map of uncertainties for nu mu events
void toyMC::makeFVMapNuE(int nmcmcpts, const char* outfile){

  // make array of histos
  cout<<"Initializing array of histograms..."<<endl;
  TH1D* hE = new TH1D("hE_nuE_v2","hE_nuE",EnuNBinsElectron,EnuBinningElectron);
  TH2FV* hfv = new TH2FV("hfv",1);
  // array of nu energy histograms
  hArrFV = new modHistoArrayFV(hE,hfv,nmcmcpts);

  // get list of mc events
  int nevmax = chMC->GetEntries();
  if (nMCevents>nevmax) nMCevents = nevmax;

  // get list of points in mcmc parameter space
  cout<<"Making list of MCMC points"<<endl;
  randomList* mcmclist = new randomList(nmcmcpts,chPars->GetEntries(),nmcmcpts);

  // loop over mcmc points
  for (int i=0; i<nmcmcpts; i++){

    // read in parameters
    cout<<"getting event"<<mcmclist->getAt(i)<<endl;
    chPars->GetEntry(mcmclist->getAt(i));

    // modify attributes using thes parameters
    if (i>0) modifier->setFromMCMC();

    // loop over T2K MC events
    for (int iev=0; iev<nMCevents; iev++){
     
      // apply parameters and see if it passes cuts
//      int ipass = applyCutsToModifiedEvent(iev);
      // apply parameters and see if it passes cuts
      int ipass = 0;
      float ww = 1.0;
      if (i!=0){
        ipass = applyCutsToModifiedEvent(iev,true);
//        ww *= modifier->getEvtWeight( fastevents->vbin[iev], fastevents->vsample[iev], fastevents->vmode[iev], fastevents->vpmomv[iev] );
      }
      else{
        ipass = applyCutsToModifiedEvent(iev,false);
      }

      // if it passes fill histos
      if (ipass!=1) continue;
 
      // modified neutrino energy
      float enu = fastevents->vfqenue[iev];

      // fill total nev
      int fvbin = hArrFV->hFV[i]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],ww*fastevents->vweight[iev]) - 1;

      // fill FV histos for different event catagories
      int catagory = getEventCatagory(iev,12);
      //
      if (catagory==1){ hArrFV->hFVCCQE[i]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],ww*fastevents->vweight[iev]);}
      if (catagory==2){ hArrFV->hFVCCnQE[i]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],ww*fastevents->vweight[iev]);}
      if (catagory==3){ hArrFV->hFVCCWrong[i]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],ww*fastevents->vweight[iev]);}
      if (catagory==4){ hArrFV->hFVNC[i]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],ww*fastevents->vweight[iev]);}
      if (catagory==0){cout<<"event # "<<iev<<" has not a catagory..."<<endl;}


      if (fvbin>=0) hArrFV->getHistogram(i,fvbin)->Fill(enu, ww*fastevents->vweight[iev]);

    }
  }

  // calculate summary and save output
  hArrFV->calcSummary();
  hArrFV->saveSummary(outfile);
  hArrFV->saveClose();

  return;

}
*/


////////////////////////////////////////////////////////////////
// Get event catagory. Current catagories:
//   1 -> CCQE
//   2 -> CCnQE
//   3 -> CCMisID
//   4 -> NC
//   0 -> Uncatagorized
int toyMC::getEventCatagory(int iev, int inutype){
      

      // is NC?
      if (fastevents->vmode[iev]>=30){
        return 4;
      }

      //  CC?
      
      // mid-IDed
      if (fastevents->vnutype[iev]!=inutype){
        return 3;
      }
      
      // ccqe
      else if (fastevents->vmode[iev]==1) {
        return 1;
      }

      else {
        return 2;
      }
      
      //
      return 0;
}



////////////////////////////////////////////////////////////////
// make map of uncertainties in different (wall,towall) regions
void toyMC::makeFVUncMap(int nmcmcpts, int nselection, const char* outfile, int fvbintype){

  // get name of selection
  TString selection_name;
  int selection_nutype;
  if (nselection==1){
    selection_name = Form("NuE_Bintype%d",fvbintype);
    selection_nutype = 12;
  }
  else if (nselection==3){
    selection_name = Form("NuE_1Pi_Bintype%d",fvbintype);
    selection_nutype = 12;
  }
  else if (nselection==2){
    selection_name = Form("NuMu_Bintype%d",fvbintype);
    selection_nutype = 14;
  }

  // make array of histos
  cout<<"makeFVUncMap: Initializing array of histograms..."<<endl;
  TH1D* hE; //< nu energy binning...this is passed on optimusPrime
  if (nselection==2){
     hE = new TH1D(selection_name.Data(),selection_name.Data(),EnuNBins,EnuBinning);
  }
  else if (nselection==1 || nselection==3){
     hE = new TH1D(selection_name.Data(),selection_name.Data(),EnuNBinsElectron,EnuBinningElectron);
  }

  // make FV histogram for finding bins
  // can use different bin type flags
  cout<<"makeFVUncMap: Initializing seed for FV histogram arrays..."<<endl;
  TH2FV* hfv = new TH2FV(Form("hfv_bintype%d",fvbintype),fvbintype);
  hfv->Draw();

  // array of nu energy and FV histograms to be filled
  cout<<"makeFVUncMap: Initializing histogram arrays"<<endl;
  modHistoArrayFV* hArrFV = new modHistoArrayFV(hE,hfv,nmcmcpts);

  // get list of mc events
  int nevmax = chMC->GetEntries();
  if (nMCevents>nevmax) nMCevents = nevmax;

  // get list of points in mcmc parameter space
  cout<<"Making list of MCMC points"<<endl;
  randomList* mcmclist = new randomList(nmcmcpts,chPars->GetEntries(),nmcmcpts);

  // loop over mcmc points
  for (int i=0; i<nmcmcpts; i++){

    // read in parameters
    cout<<"getting event"<<mcmclist->getAt(i)<<endl;
    chPars->GetEntry(mcmclist->getAt(i));

    // modify attributes using thes parameters
    modifier->setFromMCMC();

    // loop over T2K MC events
    for (int iev=0; iev<nMCevents; iev++){

      // apply parameters and see if it passes cuts
      int ipass = 0;
      float ww = 1.0;
      if (i!=0){
        ipass = applyCutsToModifiedEvent(iev,true); 
        ww *= modifier->getEvtWeight( fastevents->vbin[iev], fastevents->vsample[iev], fastevents->vmode[iev], fastevents->vpmomv[iev] );
      }
      else{
        ipass = applyCutsToModifiedEvent(iev,false);
      }

      // if it passes fill histos
      if (ipass!=nselection) continue;

      // modified nu energy
      float enu = fastevents->vfqenumu[iev];

      // fill total nev
      int fvbin = hArrFV->hFV[i]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],ww*fastevents->vweight[iev]) - 1;

      // fill FV histos for different event catagories
      //   1 -> CCQE
      //   2 -> CCnQE
      //   3 -> CCMisID
      //   4 -> NC
      //   0 -> Uncatagorized
      int catagory = getEventCatagory(iev,selection_nutype);
      // fill histograms based on catagory
      if (catagory==1){ hArrFV->hFVCCQE[i]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],ww*fastevents->vweight[iev]);}
      if (catagory==2){ hArrFV->hFVCCnQE[i]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],ww*fastevents->vweight[iev]);}
      if (catagory==3){ hArrFV->hFVCCWrong[i]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],ww*fastevents->vweight[iev]);}
      if (catagory==4){ hArrFV->hFVNC[i]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],ww*fastevents->vweight[iev]);}
      if (catagory==0){cout<<"event # "<<iev<<" has not a catagory..."<<endl;}
      if (fvbin>=0) hArrFV->getHistogram(i,fvbin)->Fill(enu, ww*fastevents->vweight[iev]);

    }
  }

  // calculate summary and save output
  hArrFV->calcSummary();
  hArrFV->saveSummary(outfile);
  hArrFV->saveClose();
  return;



}

/*
////////////////////////////////////////////////////////////////
// make map of uncertainties for nu mu events
void toyMC::makeFVMapNuMu(int nmcmcpts, const char* outfile){

  // make array of histos
  cout<<"Initializing array of histograms..."<<endl;
  TH1D* hE = new TH1D("hE_nuMu","hE_nuMu",EnuNBins,EnuBinning);
//  TH1D* hE = new TH1D("hE_nuMu","hE_nuMu",20,0,3000);

  TH2FV* hfv = new TH2FV("hfv",2);
  // array of nu energy histograms
  hArrFV = new modHistoArrayFV(hE,hfv,nmcmcpts);

  // get list of mc events
  int nevmax = chMC->GetEntries();
  if (nMCevents>nevmax) nMCevents = nevmax;

  // get list of points in mcmc parameter space
  cout<<"Making list of MCMC points"<<endl;
  randomList* mcmclist = new randomList(nmcmcpts,chPars->GetEntries(),nmcmcpts);

  // loop over mcmc points
  for (int i=0; i<nmcmcpts; i++){

    // read in parameters
    cout<<"getting event"<<mcmclist->getAt(i)<<endl;
    chPars->GetEntry(mcmclist->getAt(i));

    // modify attributes using thes parameters
    modifier->setFromMCMC();

    // loop over T2K MC events
    for (int iev=0; iev<nMCevents; iev++){

      // apply parameters and see if it passes cuts
      int ipass = 0;
      float ww = 1.0;
      if (i!=0){
        ipass = applyCutsToModifiedEvent(iev,true); 
//        ww *= modifier->getEvtWeight( fastevents->vbin[iev], fastevents->vsample[iev], fastevents->vmode[iev], fastevents->vpmomv[iev] );
      }
      else{
        ipass = applyCutsToModifiedEvent(iev,false);
      }

      // if it passes fill histos
      if (ipass!=2) continue;

      // modified nu energy
      float enu = fastevents->vfqenumu[iev];

      // fill total nev
      int fvbin = hArrFV->hFV[i]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],ww*fastevents->vweight[iev]) - 1;

      // fill FV histos for different event catagories
      int catagory = getEventCatagory(iev,14);

      if (catagory==1){ hArrFV->hFVCCQE[i]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],ww*fastevents->vweight[iev]);}
      if (catagory==2){ hArrFV->hFVCCnQE[i]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],ww*fastevents->vweight[iev]);}
      if (catagory==3){ hArrFV->hFVCCWrong[i]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],ww*fastevents->vweight[iev]);}
      if (catagory==4){ hArrFV->hFVNC[i]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],ww*fastevents->vweight[iev]);}
      if (catagory==0){cout<<"event # "<<iev<<" has not a catagory..."<<endl;}


      if (fvbin>=0) hArrFV->getHistogram(i,fvbin)->Fill(enu, ww*fastevents->vweight[iev]);
    }
  }

  // calculate summary and save output
  hArrFV->calcSummary();
  hArrFV->saveSummary(outfile);
  hArrFV->saveClose();

  return;

}
*/

/*
/////////////////////////////////////////////////////////////////
// see uncertainty in reconstructed enu spectrum
void toyMC::runToyNuMuEnu(int nmcmcpts, int nmcevents){

   // make array of histos
   cout<<"Initializing array of histograms..."<<endl;
   TH1D* hE = new TH1D("hE","hE",20,0,2000);
   hArr = new modHistoArray(hE,nmcmcpts);

   // get list of mc events
   int nevmax = chMC->GetEntries();
   cout<<"Making list of MC events"<<endl;
   randomList* evlist = new randomList(nmcevents,nevmax,nmcevents);
  
   // for selecting signal
   eventSelector* evsel = new eventSelector();
  
   // get list of mcmc points
   cout<<"Making list of MCMC points"<<endl;
   randomList* mcmclist = new randomList(nmcmcpts,chPars->GetEntries(),nmcmcpts);

   // loops
   for (int i=0; i<nmcmcpts; i++){

     // read MCMC pars
     cout<<"point: "<<mcmclist->getAt(i)<<endl;
     chPars->GetEntry(mcmclist->getAt(i)); //< read in post-fit pars
     modifier->setFromMCMC(); //< set parameters 

     for (int iev=0; iev<nmcevents; iev++){
       
       // read in event
       if ((iev%100)==0){
         cout<<"getting entry "<<iev<<endl;
         cout<<"getting entry "<<evlist->getAt(iev)<<endl;
       }

//       chMC->GetEntry(evlist->getAt(iev));
       chMC->GetEntry((iev));
       if (mcEvent->ncomponent==5) continue;
       modifier->applyPars(); //< actually modifies att[]

       // see if event passes cuts
       double lmom = mcEvent->attribute[3];
       double lpid = mcEvent->attribute[0];
       int ipass = evsel->selectNuMu(mcEvent->nhitac,
                                     lmom,
                                     lpid,
                                     0,0,0);
       // if it passes fill histos
       if (ipass!=1) continue;
       TVector3 rcdir;
       rcdir.SetXYZ(mcEvent->fq1rdir[0][2][0],mcEvent->fq1rdir[0][2][1],mcEvent->fq1rdir[0][2][2]);
       float enu_rc = calcEnu(lmom,rcdir,1);
       hArr->histos[i]->Fill(enu_rc,mcEvent->evtweight);
     }
   }

   return;
}
*/

/*
////////////////////////////////////////////////////////////////
// 
void toyMC::testToy(int nmcmcpts){
 

   // test seed histograms
//   TH1D* hE = new TH1D("hE","hE",20,0,2000);
//   hArr = new modHistoArray(hE,nmcmcpts);
   TH1D* hPID = new TH1D("hpid","hpid",30,-2000,2000);
   hArr = new modHistoArray(hPID,nmcmcpts);

   // histo array

   // max # of events
   int nevmax = chMC->GetEntries();

   // for selecting signal
   eventSelector* evsel = new eventSelector();
  
   // get list of mcmc points
   int thepoints[100];
   for (int i=0; i<nmcmcpts; i++){
     thepoints[i] =  getRandomMCMCPoint();
   }
  
   // loops
   nevmax = 10000; //< max number of mc events to use
   for (int i=0; i<nmcmcpts; i++){
     cout<<"point: "<<i<<endl;
     chPars->GetEntry(thepoints[i]); //< read in post-fit pars
     modifier->setFromMCMC(); //< set parameters 
     for (int iev=0; iev<nevmax; iev++){
       if ((iev%100)==0) cout<<"getting entry "<<iev<<endl;
       chMC->GetEntry(iev);
       modifier->applyPars(); //< actually modifies att[]

       // see if event passes cuts
       double lmom = mcEvent->attribute[3];
       double lpid = mcEvent->attribute[0];
//       int ipass = evsel->selectNuMu(mcEvent->nhitac,
//                                     lmom,
//                                     lpid,
//                                     0,0,0);
//       if (ipass!=1) continue;
//
       hArr->histos[i]->Fill(mcEvent->attribute[0]);
     }
   }


   return;
}
*/




void toyMC::setAtmFitPars(const char* parfile){
  
  fitPars = new atmFitPars(parfile); 

  modifier = new mcmcApply(fitPars, mcmcPars);

}


////////////////////////////////////////////////////////////////
void toyMC::setCompare(histoCompare* hc){

 hCompare = hc;

 modifier = new mcmcApply(hCompare->thePars, mcmcPars);

 return;
}



////////////////////////////////////////////////////////////////
/*
void toyMC::fillArrayDirect(int isamp, int ibin, int iatt, int npts){

  // get seed
  TH1D* hseed = hCompare->hManager->getSumHistogramMod(isamp,ibin,iatt);
  hseed->Draw();

  hArr = new modHistoArray(hseed,npts);

  int nmcmcpts = chPars->GetEntries();
  for (int ipt=0; ipt<npts; ipt++){
    int randpt = randy->Integer(nmcmcpts);
    mcmcPars->GetEntry(randpt);
    modifier->setFromMCMC();
    TH1D* hadd = hCompare->hManager->getSumHistogramMod(isamp,ibin,iatt);
    hArr->setHistoContents(hadd);
  }

  hCompare->hManager->hData[isamp][ibin][iatt]->SetMarkerStyle(8);
  hCompare->hManager->hData[isamp][ibin][iatt]->SetMarkerSize(1.2);
  hCompare->hManager->hData[isamp][ibin][iatt]->Draw("e");
  hCompare->hManager->getSumHistogram(isamp,ibin,iatt)->SetLineColor(kRed);
  hCompare->hManager->getSumHistogram(isamp,ibin,iatt)->SetLineWidth(2);
  hCompare->hManager->getSumHistogram(isamp,ibin,iatt)->Draw("sameh");
  for (int ipt=0; ipt<npts; ipt++){
    hArr->histos[ipt]->SetLineColor(kBlue);
    hArr->histos[ipt]->SetLineWidth(1);
    hArr->histos[ipt]->Draw("sameh");
  }

  return;
}

*/
















#endif
