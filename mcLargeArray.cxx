#ifndef LARGEARRAY_CXX
#define LARGEARRAY_CXX


#include "mcLargeArray.h"

mcLargeArray::mcLargeArray(TChain* ch, int nevents){
 
  nsize = nevents;
  fillArray(ch);
//  fillThinArray(ch,2);
}

void mcLargeArray::fillThinArray(TChain* ch, int thinning){


 fqProcessedEvent* mcevent = new fqProcessedEvent(ch);

 int nfilled = 0;
 for (int i=0; i<ch->GetEntries(); i++){
  int ievent = i;
  if ((ievent%thinning)==0) continue;
  ch->GetEntry(ievent);
  nfilled++;
  vnutype[i] = (Short_t)mcevent->ipnu[0];
  vfqmumom[i] = (float)mcevent->fq1rmom[0][2];
  vfqemom[i] = (float)mcevent->fq1rmom[0][1];
  vnhitac[i] = mcevent->nhitac;
  vfqpi0par[i] = (float)mcevent->fqpi0par;
  vfqwall[i] = (float)mcevent->fqwall;
  vfqtowall[i] = (float)mcevent->fqtowall;
  vwallv[i] = (float)mcevent->wallv;
  vmode[i] = (Short_t)TMath::Abs(mcevent->mode);
  voscpower[i][0] = (float)mcevent->oscpower[0];
  voscpower[i][1] = (float)mcevent->oscpower[1];
  voscpower[i][2] = (float)mcevent->oscpower[2];
  voscpower[i][3] = (float)mcevent->oscpower[3];
  vfqnring[i] = (int)mcevent->fqmrnring[0];
  vfqpid[i] = (float)mcevent->fq1rnll[0][2]-(float)mcevent->fq1rnll[0][1];
  vweight[i] = (float)mcevent->evtweight;
  vfqenue[i] = (float)mcevent->fq1renu[0];
  vfqenumu[i] = (float)mcevent->fq1renu[1];
  vfqnsubev[i] = (int)mcevent->fqnse;
 }

 nsize = nfilled;

 return;
}


void mcLargeArray::fillArray(TChain* ch){


 fqProcessedEvent* mcevent = new fqProcessedEvent(ch);

 TRandom2* randy = new TRandom2(nsize);
 int nmax = ch->GetEntries();
 int flgUseRandom=1;
 if (nsize>=nmax){ 
   flgUseRandom=0;
   nsize=nmax;
 }

//  fill array with random events
// const int NN = nsize;
// int Events[NN];
// for (int i=0; i<nsize; i++){
//   Events[i] = randy->Integer(nmax);
// }
// sort(Events,Events+NN);

 for (int i=0; i<nsize; i++){
  int ievent = i;
//  if (flgUseRandom) ievent = Events[i]; 
  if ((i%100)==0) cout<<ievent<<endl;
  ch->GetEntry(ievent);
  vnutype[i] = (Short_t)mcevent->ipnu[0];
  vfqmumom[i] = (float)mcevent->fq1rmom[0][2];
  vfqemom[i] = (float)mcevent->fq1rmom[0][1];
  vnhitac[i] = mcevent->nhitac;
  vfqpi0par[i] = (float)mcevent->fqpi0par;
  vfqwall[i] = (float)mcevent->fqwall;
  vfqtowall[i] = (float)mcevent->fqtowall;
  vwallv[i] = (float)mcevent->wallv;
  vmode[i] = (Short_t)mcevent->mode;
  voscpower[i][0] = (float)mcevent->oscpower[0];
  voscpower[i][1] = (float)mcevent->oscpower[1];
  voscpower[i][2] = (float)mcevent->oscpower[2];
  voscpower[i][3] = (float)mcevent->oscpower[3];
  vfqnring[i] = (int)mcevent->fqmrnring[0];
  vfqpid[i] = (float)mcevent->fq1rnll[0][2]-(float)mcevent->fq1rnll[0][1];
  vweight[i] = (float)mcevent->evtweight;
 }

 return;
}

#endif
