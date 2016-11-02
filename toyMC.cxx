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
// 
void toyMC::testToy(int nmcmcpts){
 

   // test seed histograms
   TH1D* hE = new TH1D("hE","hE",20,0,2000);
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
   nevmax = 1000; //< max number of mc events to use
   for (int i=0; i<nmcmcpts; i++){
     cout<<"point: "<<i<<endl;
     chPars->GetEntry(thepoints[i]); //< read in post-fit pars
     modifier->setFromMCMC(); //< set parameters 
     for (int iev=0; iev<nevmax; iev++){
       if ((iev%100)==0) cout<<"getting entry "<<iev<<endl;
       chMC->GetEntry(iev);
       modifier->applyPars(); //< actually modifies att[]
       hArr->histos[i]->Fill(mcEvent->attribute[0]);
     }
   }


   return;
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
