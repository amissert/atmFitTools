#ifndef PRIME_CXX
#define PRIME_CXX


#include "optimusPrime.h"

//////////////////////////////////////////
// constructor
optimusPrime::optimusPrime(TChain* t2kmc, int nevts){

 chmc = t2kmc;
 mcevent = new fqProcessedEvent(chmc);
 nevents = nevts;
 FOMType = 0;

 if (nevts>chmc->GetEntries()){
   nevents = chmc->GetEntries();
   flgUseEventList = 1;
 }
 else{
  flgUseEventList = 1;
 }

 // we only need enu, wall, towall, and mode and wallv
 chmc->SetBranchStatus("*",0);
 chmc->SetBranchStatus("mode",1);
 chmc->SetBranchStatus("wallv",1);
 chmc->SetBranchStatus("nhitac",1);
 chmc->SetBranchStatus("attribute",1);
 chmc->SetBranchStatus("oscpower",1);
 chmc->SetBranchStatus("evtweight",1);
 chmc->SetBranchStatus("fq*",1);
 chmc->SetBranchStatus("ipnu",1);

 hFVAll = new TH2FV("hall",-1,30,0,800,30,0,800);
 hFVAvg = new TH2FV("havg",-1,30,0,800,30,0,800);

 uncertaintyCalculator = new moreUncertainties("/nfs/data41/t2k/amissert/atmos/head/atmFitTools/data/");

 fillArray();

}

/////////////////////////////////////////
// get any additional uncertainty
float optimusPrime::getMoreUncertainty(float wallv, float wallrc){
  return uncertaintyCalculator->getTotalUncertainty(wallv,wallrc);
}

/////////////////////////////////////////
// get the oscillation power for this event
float optimusPrime::getOscPower(int nutype, int oscpar){

 // NC events do not contribute 
 if (TMath::Abs(mcevent->mode)>=30) return 0.;
 
 // other nu do not contribute
 if (TMath::Abs(mcevent->ipnu[0]!=nutype)) return 0.;

 float oscpow = (float)mcevent->evtweight*(float)mcevent->oscpower[oscpar]; 

 return oscpow;

}

/////////////////////////////////////////
// get the oscillation power for this event
float optimusPrime::getOscPowerFast(int nutype, int ientry, int oscpar){

 // NC events do not contribute 
 if (TMath::Abs(fastevents->vmode[ientry])>=30) return 0.;
 
 // other nu do not contribute
 if (TMath::Abs(fastevents->vnutype[ientry] != nutype)) return 0.;

 float oscpow = fastevents->vweight[ientry]*fastevents->voscpower[ientry][oscpar]; 

 return oscpow;

}

void optimusPrime::fillFVHistoFast(){
//  hFVAll = new TH2FV("hall",-1,40,0,1000,40,0,1000);
//  hFVAvg = new TH2FV("havg",-1,40,0,1000,40,0,1000);
  hFVAll->Reset();
  hFVAvg->Reset();
  for (int ievt=0; ievt<nevents; ievt++){
    hFVAll->Fill(fastevents->vfqtowall[ievt],fastevents->vfqwall[ievt],fastevents->vweight[ievt]);
    if (AvgType==0){
       // Signal
       if ((fastevents->vmode[ievt]==1)&&(fastevents->vnutype[ievt]==14)){
         hFVAvg->Fill(fastevents->vfqtowall[ievt],fastevents->vfqwall[ievt],fastevents->vweight[ievt]);
       }
    }
    if (AvgType==1){
       // BG 
       if (!((fastevents->vmode[ievt]==1)&&(fastevents->vnutype[ievt]==14))){
         hFVAvg->Fill(fastevents->vfqtowall[ievt],fastevents->vfqwall[ievt],fastevents->vweight[ievt]);
       } 
    }
    if (AvgType==2){    
      // Power 
      hFVAvg->Fill(fastevents->vfqtowall[ievt],fastevents->vfqwall[ievt],fastevents->vweight[ievt]*fastevents->voscpower[ievt][0]);      
    }
  }
  hFVAll->Draw("colz");
}

void optimusPrime::calcFOMMap(float towallmax, float wallmax,int oscpar){
 
 hFV = new TH2FV("h",-1,30,0,towallmax,30,0,wallmax);
 hFV->SetContour(100);
 for (int ibin = 0; ibin<hFV->GetNumberOfBins(); ibin++){
   float wallcut = (float)hFV->GetBinCenterY(ibin);
   float towallcut = (float)hFV->GetBinCenterX(ibin);
   float value = calcNuMuFOM(towallcut,wallcut,oscpar);
   hFV->SetBinContent(ibin,value);
 }
 hFV->Draw("colz");
  return;
}

//////////////////////////////////////////
//Read events into memory for fast looping
void optimusPrime::fillArray(){

  fastevents = new mcLargeArray(chmc,nevents); 

  //
  return;
}


/////////////////////////////////////////
// get total FOM
float optimusPrime::calcNuMuFOM(float towallmin, float wallmin, int oscpar){

  if (towallmin<wallmin) return 0.;

  Nevents = 0.;
  NB = 0.; 
  NS = 0.; 
  Power = 0.;
  Syst = 0.;

  // use short event list
//  if (flgUseEventList){
    for (int i=0; i<nevents; i++){
//      if ((i%50000)==0) cout<<i<<endl;
      int ipass = selectNuMu(fastevents->vnhitac[i],
                             fastevents->vfqmumom[i],
                             fastevents->vfqpid[i],
                             fastevents->vfqnring[i]);
       if (ipass){
        if ((fastevents->vfqwall[i] > wallmin)&&(fastevents->vfqtowall[i]>towallmin)){
          Power += getOscPowerFast(14,i,oscpar);
          Nevents += fastevents->vweight[i];
//          Syst += getMoreUncertainty(fastevents->vwallv[i],fastevents->vfqwall[i])*fastevents->vweight[i];
          Syst += uncertaintyCalculator->getTotalUncertainty(fastevents->vwallv[i],fastevents->vfqwall[i])*fastevents->vweight[i];
          if ((fastevents->vmode[i]!=1 || fastevents->vnutype[i]!=14)){
            NB+=fastevents->vweight[i];
          }
          else{
            NS+=fastevents->vweight[i];
          }
        }
      }
    }
//  }

  // use all events
//  else{
//    for (int i=0; i<chmc->GetEntries(); i++){
//      if ((i%1000)==0) cout<<i<<endl;
//      chmc->GetEntry(i);
//      if (mcevent->passMuCuts()){
//        if ((mcevent->fq1rwall[0][2] > wallmin)&&(mcevent->fq1rtowall[0][2]>towallmin)){
//          Power += getOscPower(14,oscpar);
//          Nevents += mcevent->evtweight;
//          Syst += 0.;
//        }
//      }
//    }
//  }
 
  float fom = (Power*Power)/(Nevents+(Syst*Syst));
  if (FOMType==1) fom = Nevents;
  if (FOMType==2) fom = Syst*Syst;
  if (FOMType==3) fom = Power*Power/(Nevents*Nevents);
  if (FOMType==4) fom = Power*Power;
  if (FOMType==5) fom = NS;
  if (FOMType==6) fom = NB;
  if (FOMType==7) fom = NS/NB;
  if (FOMType==8) fom = NS/TMath::Sqrt(NB+NS);
  if (FOMType==9) fom = NS/TMath::Sqrt(Nevents);



  float scale = (float)nevents/(float)chmc->GetEntries();
  fom*=1./scale;
  cout<<"FOM: "<<fom<<endl;
  cout<<"N: "<<Nevents<<endl;
  cout<<"S: "<<Syst<<endl;
  cout<<"P: "<<Power<<endl;

  return fom;

}

#endif







