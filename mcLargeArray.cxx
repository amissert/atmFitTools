#ifndef LARGEARRAY_CXX
#define LARGEARRAY_CXX


#include "mcLargeArray.h"

mcLargeArray::mcLargeArray(TChain* ch, int nevents){
 
  // select random events if not using everything
  nsize = nevents;
  int nmax = ch->GetEntries();
  if (nsize>=nmax){ 
    nsize=nmax;
  }
  fillArray(ch);
}

/*
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
  vbin[i] = (Short_t)mcevent->nbin;
  vcomponent[i] = (Short_t)mcevent->ncomponent;
  vattribute[i][0] = (float)mcevent->attribute[0];
  vattribute[i][1] = (float)mcevent->attribute[1];
  vattribute[i][2] = (float)mcevent->attribute[2];
  vattribute[i][3] = (float)mcevent->attribute[3];
  vfqdir[i][0][0] = (float)mcevent->fq1rdir[0][1][0];
  vfqdir[i][0][1] = (float)mcevent->fq1rdir[0][1][1];
  vfqdir[i][0][2] = (float)mcevent->fq1rdir[0][1][2];
  vfqdir[i][1][0] = (float)mcevent->fq1rdir[0][2][0];
  vfqdir[i][1][1] = (float)mcevent->fq1rdir[0][2][1];
  vfqdir[i][1][2] = (float)mcevent->fq1rdir[0][2][2];
 }

 nsize = nfilled;

 return;
}
*/

void mcLargeArray::fillArray(TChain* ch){

 ch->SetBranchStatus("*",0);
 ch->SetBranchStatus("fq*",1);
 ch->SetBranchStatus("*v",1);
 ch->SetBranchStatus("nhitac",1);
 ch->SetBranchStatus("evtweight",1);
 ch->SetBranchStatus("oscpower",1);
 ch->SetBranchStatus("mode",1);
 ch->SetBranchStatus("ipnu",1);
 ch->SetBranchStatus("pass*",1);
 ch->SetBranchStatus("nbin",1);
 ch->SetBranchStatus("nsample",1);
 ch->SetBranchStatus("ncomponent",1);
 ch->SetBranchStatus("att*",1);

 

 fqProcessedEvent* mcevent = new fqProcessedEvent(ch);



 for (int i=0; i<nsize; i++){
  int ievent = i;
//  if (flgUseRandom) ievent = Events[i]; 
  if ((i%10000)==0) cout<<"filling large array...on event: "<<ievent<<endl;
  ch->GetEntry(ievent);
  vnutype[i] = TMath::Abs((Short_t)mcevent->ipnu[0]);
  vfqmumom[i] = (float)mcevent->fq1rmom[0][2];
  vfqemom[i] = (float)mcevent->fq1rmom[0][1];
  vnhitac[i] = mcevent->nhitac;
  vfqpi0par[i] = (float)mcevent->fqpi0par;
  vfqwall[i] = (float)mcevent->fqwall;
  vfqtowall[i] = (float)mcevent->fqtowall;
  vwallv[i] = (float)mcevent->wallv;
  vtowallv[i] = (float)mcevent->towallv[2];
  vmode[i] = (Short_t)TMath::Abs(mcevent->mode);
  voscpower[i][0] = (float)mcevent->oscpower[0];
  voscpower[i][1] = (float)mcevent->oscpower[1];
  voscpower[i][2] = (float)mcevent->oscpower[2];
  voscpower[i][3] = (float)mcevent->oscpower[3];
  vfqnring[i] = (int)mcevent->fqmrnring[0];
  vfqpid[i] = (float)mcevent->fq1rnll[0][2]-(float)mcevent->fq1rnll[0][1];
  vweight[i] = (float)mcevent->evtweight; 
  vbin[i] = (Short_t)mcevent->nbin;
  vsample[i] = (Short_t)mcevent->nsample;
  vfqnsubev[i] = (int)mcevent->fqnse;
  vcomponent[i] = (Short_t)mcevent->ncomponent;
  vattribute[i][0] = (float)mcevent->attribute[0];
  vattribute[i][1] = (float)mcevent->attribute[1];
  vattribute[i][2] = (float)mcevent->attribute[2];
  vattribute[i][3] = (float)mcevent->attribute[3];
  vfqdir[i][0][0] = (float)mcevent->fq1rdir[0][1][0];
  vfqdir[i][0][1] = (float)mcevent->fq1rdir[0][1][1];
  vfqdir[i][0][2] = (float)mcevent->fq1rdir[0][1][2];
  vfqdir[i][1][0] = (float)mcevent->fq1rdir[0][2][0];
  vfqdir[i][1][1] = (float)mcevent->fq1rdir[0][2][1];
  vfqdir[i][1][2] = (float)mcevent->fq1rdir[0][2][2];
  vfqenue[i] = (float)mcevent->fq1renu[0];
  vfqenumu[i] = (float)mcevent->fq1renu[1];
  vfqrcpar[i] = (float)mcevent->fqrcpar;
  vpmomv[i] = (float)mcevent->pmomv[0];
  vpassnumu[i] = (Short_t)mcevent->passmucut;
  vpassnue[i] = (Short_t)mcevent->passecut;
  vpassnue1rpi[i] = (Short_t)mcevent->passe1rpicut;
 }
 
 ch->SetBranchStatus("*",1);

 return;
}

//////////////////////////////////////////////
// print event info
void mcLargeArray::printEvent(int iev){

  cout<<"-----Event # "<<iev<<"-----"<<endl;
  cout<<" bin:         "<<vbin[iev]<<endl;
  cout<<" component:   "<<vcomponent[iev]<<endl;
  cout<<" mu momentum: "<<vfqmumom[iev]<<endl;
  cout<<" e momentum:  "<<vfqemom[iev]<<endl;
  cout<<" nu momentum:  "<<vpmomv[iev]<<endl;
  cout<<" PID:         "<<vfqpid[iev]<<endl;
  cout<<" Pi0:         "<<vfqpi0par[iev]<<endl;
  cout<<" PiP:         "<<vattribute[iev][2]<<endl;
  cout<<" mode:        "<<vmode[iev]<<endl;
  cout<<" nring:        "<<vfqnring[iev]<<endl;

  return;

}
#endif
