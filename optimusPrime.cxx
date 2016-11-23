#ifndef PRIME_CXX
#define PRIME_CXX


#include "optimusPrime.h"

//////////////////////////////////////////
// constructor
optimusPrime::optimusPrime(TChain* t2kmc, int nevts){

 chmc = t2kmc;
 mcevent = new fqProcessedEvent(chmc);

 if (nevts>chmc->GetEntries()){
   nevts = chmc->GetEntries();
   flgUseEventList = 0;
 }
 else{
  flgUseEventList = 1;
  eventlist = new randomList(nevts, chmc->GetEntries(),nevts);
 }

 // we only need enu, wall, towall, and mode and wallv
 chmc->SetBranchStatus("*",0);
 chmc->SetBranchStatus("fqwall",1);
 chmc->SetBranchStatus("fqtowall",1);
 chmc->SetBranchStatus("fqenu",1);
 chmc->SetBranchStatus("mode",1);
 chmc->SetBranchStatus("wallv",1);
 chmc->SetBranchStatus("ipnu",1);


}

/////////////////////////////////////////
// get the oscillation power for this event
float optimusPrime::getOscPower(int nutype, int oscpar){

 // NC events do not contribute 
 if (TMath::Abs(mcevent->mode)>=30) return 0.;
 
 // other nu do not contribute
 if (TMath::Abs(mcevent->ipnu[0]!=nutype)) return 0.;

 float oscpow = mcevent->evtweight*mcevent->oscpower[oscpar]; 

 return oscpow;

}


/////////////////////////////////////////
// get total FOM
void optimusPrime::calcNuMuFOM(float towallmin, float wallmin, int oscpar){

  float P = 0.;
  float N = 0.;
  float S = 0.;

  // use short event list
  if (flgUseEventList){
    for (int i=0; i<eventlist->nMax; i++){
      chmc->GetEntry(eventlist->getAt(i));
      if ((i%1000)==0) cout<<i<<endl;
      if (mcevent->passMuCuts()){
        if ((mcevent->fq1rwall[0][2] > wallmin)&&(mcevent->fq1rtowall[0][2]>towallmin)){
          P += getOscPower(14,oscpar);
          N += mcevent->evtweight;
          S += 0.;
        }
      }
    }
  }

  // use all events
  else{
    for (int i=0; i<chmc->GetEntries(); i++){
      if ((i%1000)==0) cout<<i<<endl;
      chmc->GetEntry(i);
      if (mcevent->passMuCuts()){
        if ((mcevent->fq1rwall[0][2] > wallmin)&&(mcevent->fq1rtowall[0][2]>towallmin)){
          P += getOscPower(14,oscpar);
          N += mcevent->evtweight;
          S += 0.;
        }
      }
    }
  }
 
  float fom = (P*P)/(N+S);
  cout<<"FOM: "<<fom<<endl;
  cout<<"N: "<<fom<<endl;
  cout<<"S: "<<fom<<endl;
  cout<<"P: "<<fom<<endl;

  return;

}

#endif







