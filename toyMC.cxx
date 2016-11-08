#ifndef TOYMC_CXX
#define TOYMC_CXX

#include "toyMC.h"

using namespace std;

/////////////////////////////////////////////////////////////////
toyMC::toyMC(){


}

/////////////////////////////////////////////////////////////////
// set the pointers to the mc events and the post-mcmc parameters
void toyMC::setChains(TChain* chmc, TChain* chpars){

  chMC = chmc;
  chPars = chpars;

  mcEvent = new fqProcessedEvent(chMC);
  mcEvent->useImportantOnly();

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


////////////////////////////////////////////////////////////////
// make map of uncertainties for nu mu events
void toyMC::makeFVMapNuMu(int nmcmcpts, int mcevents){

   // make array of histos
   cout<<"Initializing array of histograms..."<<endl;
   TH1D* hE = new TH1D("hE","hE",20,0,2000);
   TH2FV* hfv = new TH2FV("hfv",1);
   hArrFV = new modHistoArrayFV(hE,hfv,nmcmcpts);

   // get list of mc events
   int nevmax = chMC->GetEntries();
   if (mcevents<nevmax) nevmax = mcevents;

   // for selecting signal
   eventSelector* evsel = new eventSelector();
  
   // get list of mcmc points
   cout<<"Making list of MCMC points"<<endl;
   randomList* mcmclist = new randomList(nmcmcpts,chPars->GetEntries(),nmcmcpts);

   /*
   // loops!
   for (int iev=0; iev<nevmax; iev++){

   // read in event
   if ((iev%100)==0) cout<<"toyMC: getting event: "<<iev<<endl;
   chMC->GetEntry(iev);
   if (mcEvent->ncomponent==5) continue;
   
   hfv->Reset();

   for (int i=0; i<nmcmcpts; i++){
     
     // read mcmc pars and apply them
     chPars->GetEntry(mcmclist->getAt(i));
     modifier->applyPars();

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
     int fvbin = hfv->Fill(mcEvent->fq1rtowall[0][2],mcEvent->fq1rwall[0][2]) - 1;
     if (fvbin>=0) hArrFV->getHistogram(i,fvbin)->Fill(enu_rc, mcEvent->evtweight);

    } //< end MCMC pars loop

  } //< end MC events loop
  */

  for (int i=0; i<nmcmcpts; i++){
    chPars->GetEntry(mcmclist->getAt(i));
    modifier->setFromMCMC();
    for (int iev=0; iev<nevmax; iev++){
      // read in event
      if ((iev%5000)==0) cout<<"toyMC: getting event: "<<iev<<endl;
      chMC->GetEntry(iev);
      if (mcEvent->ncomponent==5) continue;
      modifier->applyPars();
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
      int fvbin = hArrFV->hFV[i]->Fill(mcEvent->fq1rtowall[0][2],mcEvent->fq1rwall[0][2],mcEvent->evtweight) - 1;
      if (fvbin>=0) hArrFV->getHistogram(i,fvbin)->Fill(enu_rc, mcEvent->evtweight);
    }
  }
  return;

}



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





















#endif