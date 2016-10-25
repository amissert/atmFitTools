#ifndef TOYMC_CXX
#define TOYMC_CXX

#include "toyMC.h"

using namespace std;

/////////////////////////////////////////////////////////////////
toyMC::toyMC(){


}

/////////////////////////////////////////////////////////////////
void toyMC::setChains(TChain* chmc, TChain* chpars){

  chMC = chmc;
  chPars = chpars;

  mcEvent = new fqProcessedEvent(chMC);
  mcmcPars = new mcmcReader(chPars);

  return;
}

////////////////////////////////////////////////////////////////
int toyMC::getRandomMCMCPoint(){
  int nmax = chPars->GetEntries();
  int randpoint = randy->Integer(nmax);
  cout<<"getting mcmc point: "<<randpoint<<endl;
  chPars->GetEntry(randpoint);
  return randpoint;
}

////////////////////////////////////////////////////////////////
void toyMC::testToy(int nmcmcpts){
  
   TH1D* hE = new TH1D("hE","hE",20,0,2000);
   hArr = new modHistoArray(hE,nmcmcpts);
   int nevmax = chMC->GetEntries();
   eventSelector* evsel = new eventSelector();

   
   int thepoints[100];
   for (int i=0; i<nmcmcpts; i++){
     thepoints[i] =  getRandomMCMCPoint();
   }
   
//   nevmax = 1000;
//   for (int iev=0; iev<nevmax; iev++){
//     if ((iev%100)==0) cout<<"getting entry "<<iev<<endl;
//     chMC->GetEntry(iev);
//     for (int i=0; i<nmcmcpts; i++){
//       cout<<"point: "<<i<<endl;
//       chPars->GetEntry(thepoints[i]);
//       modifier->setFromMCMC(); //< modifies attribute[]
//       double pmom = mcEvent->fq1rmom[0][2];
//       double elike = mcEvent->attribute[0];
//       double nring = (double)mcEvent->fqmrnring[0];
//       if (evsel->selectNuMu(pmom,elike,nring)){
//         hArr->histos[i]->Fill(pmom);
//       }
//     }
//   }

   nevmax = 1000;
   for (int i=0; i<nmcmcpts; i++){
     cout<<"point: "<<i<<endl;
     chPars->GetEntry(thepoints[i]);
     modifier->setFromMCMC(); //< modifies attribute[]
     for (int iev=0; iev<nevmax; iev++){
       if ((iev%100)==0) cout<<"getting entry "<<iev<<endl;
       chMC->GetEntry(iev);
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

  hCompare->hManager->hData[isamp][ibin][iatt]->Draw("e");
  hCompare->hManager->getSumHistogram(isamp,ibin,iatt)->Draw("same");
  for (int ipt=0; ipt<npts; ipt++){
    hArr->histos[ipt]->SetLineColor(kBlue);
    hArr->histos[ipt]->Draw("same");
  }

  return;
}























#endif
