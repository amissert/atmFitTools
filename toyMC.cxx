#ifndef TOYMC_CXX
#define TOYMC_CXX

#include "toyMC.h"

using namespace std;


/////////////////////////////////////////////////////////////////
// apply the cuts to a modified event and see if it passes
int toyMC::applyCutsToModifiedEvent(int iev){

  // fill tmp array with "nominal" MC values
  const int natt = 4;
  float attributesTmp[natt];
  for (int iatt=0; iatt<natt; iatt++){
    attributesTmp[iatt] = fastevents->vattribute[iev][iatt];   
  }
 
  // modify tmp array by applying the histogram shape parameters
  modifier->applyPars(fastevents->vbin[iev],
                      fastevents->vcomponent[iev],
                      attributesTmp,
                      natt);

  // fill cut parameter structure using modified attributes
  if (indexPIDPar>=0) cutPars.fqpid = attributesTmp[indexPIDPar];
  if (indexPi0Par>=0) cutPars.fqpi0par = attributesTmp[indexPi0Par];
  if (indexPiPPar>=0) cutPars.fqpippar = attributesTmp[indexPiPPar];
  if (indexRCPar>=0) cutPars.fqrcpar = attributesTmp[indexRCPar];

  // other cut pars that are not modified
  cutPars.fqmome = fastevents->vfqmumom[iev];
  cutPars.fqmommu = fastevents->vfqemom[iev];
  cutPars.nhitac = fastevents->vnhitac[iev];
  cutPars.fqnsubev = fastevents->vfqnsubev[iev];
  cutPars.fqenue = fastevents->vfqenue[iev];
  cutPars.fqenumu = fastevents->vfqenumu[iev];

  // see if it passes cuts
  int passnue = selectNuE(cutPars);
  int passnumu = selectNuMu(cutPars);
  
  //
  if (passnue>0) return 1;
  if (passnumu>0) return 2;
  return 0;
  
}


/////////////////////////////////////////////////////////////////
toyMC::toyMC(){

//  indexPIDPar = 0;
//  indexPi0Par = 1;
//  indexPiPPar = 2;
//  indexRCPar  = 3;
//  indexMom    = 4;

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
      if (ipass==1) t2kToys->hEnuElectron->Fill(cutPars.fqenue, fastevents->vweight[iev]);
      if (ipass==2) t2kToys->hEnuMuon->Fill(cutPars.fqenumu, fastevents->vweight[iev]);
    }
    t2kToys->finishToyRun();
  }
  
  t2kToys->calcUncertainties();

  return;

}



////////////////////////////////////////////////////////////////
// make map of uncertainties for nu mu events
void toyMC::makeFVMapNuE(int nmcmcpts, const char* outfile){

  // make array of histos
  cout<<"Initializing array of histograms..."<<endl;
  TH1D* hE = new TH1D("hE_nuE","hE_nuE",EnuNBinsElectron,EnuBinningElectron);
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
      int ipass = applyCutsToModifiedEvent(iev);

      // if it passes fill histos
      if (ipass!=1) continue;
 
      // modified neutrino energy
      float enu = cutPars.fqenue;

      // fill total nev
      int fvbin = hArrFV->hFV[i]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],fastevents->vweight[iev]) - 1;

      // is NC?
      if (fastevents->vmode[iev]>=30){ hArrFV->hFVNC[i]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],fastevents->vweight[iev]);}

      //  CC?
      if (fastevents->vnutype[iev]!=12){hArrFV->hFVCCWrong[i]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],fastevents->vweight[iev]);}
      else if (fastevents->vmode[iev]==1) {
        hArrFV->hFVCCQE[i]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],fastevents->vweight[iev]);
      }
      else {
        hArrFV->hFVCCnQE[i]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],fastevents->vweight[iev]);
      }

      if (fvbin>=0) hArrFV->getHistogram(i,fvbin)->Fill(enu, fastevents->vweight[iev]);

    }
  }

  // calculate summary and save output
  hArrFV->calcSummary();
  hArrFV->saveSummary(outfile);
  hArrFV->saveClose();

  return;

}




////////////////////////////////////////////////////////////////
// make map of uncertainties for nu mu events
void toyMC::makeFVMapNuMu(int nmcmcpts, const char* outfile){

  // make array of histos
  cout<<"Initializing array of histograms..."<<endl;
  TH1D* hE = new TH1D("hE_nuMu","hE_nuMu",EnuNBins,EnuBinning);
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
      int ipass = applyCutsToModifiedEvent(iev);
         
      // if it passes fill histos
      if (ipass!=2) continue;

      // modified nu energy
      float enu = cutPars.fqenumu;

      // fill total nev
      int fvbin = hArrFV->hFV[i]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],fastevents->vweight[iev]) - 1;

      // is NC?
      if (fastevents->vmode[iev]>=30){ hArrFV->hFVNC[i]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],fastevents->vweight[iev]);}

      //  CC?
      if (fastevents->vnutype[iev]!=14){hArrFV->hFVCCWrong[i]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],fastevents->vweight[iev]);}
      else if (fastevents->vmode[iev]==1) {
        hArrFV->hFVCCQE[i]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],fastevents->vweight[iev]);
      }
      else {
        hArrFV->hFVCCnQE[i]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],fastevents->vweight[iev]);
      }

      if (fvbin>=0) hArrFV->getHistogram(i,fvbin)->Fill(enu, fastevents->vweight[iev]);
    }
  }

  // calculate summary and save output
  hArrFV->calcSummary();
  hArrFV->saveSummary(outfile);
  hArrFV->saveClose();

  return;

}


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

  modifier = new mcmcApply(fitPars, mcmcPars, mcEvent);

}

////////////////////////////////////////////////////////////////
void toyMC::setCompare(histoCompare* hc){

 hCompare = hc;

 modifier = new mcmcApply(hCompare->thePars, mcmcPars, mcEvent);

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
