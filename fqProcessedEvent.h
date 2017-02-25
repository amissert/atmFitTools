//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Apr 29 14:46:21 2016 by ROOT version 5.28/00c
// from TTree h1/DST
// found on file: t2kmc_numu_ppmc_2_.root
//////////////////////////////////////////////////////////

#ifndef fqProcessedEvent_h
#define fqProcessedEvent_h

#include <TROOT.h>
#include <iostream>
#include <TChain.h>
#include <TFile.h>

class fqProcessedEvent {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // index of important attribtes
   int momentumIndex;
   int PIDIndex;
   
   // Declaration of leaf types
   Int_t           nev;
   UInt_t          nhitac;
   Int_t           npar;
   Float_t         wallv;
   UInt_t          ipv[50];   //[npar]
   Float_t         posv[3];
   Float_t         dirv[50][3];   //[npar]
   Float_t         pmomv[50];   //[npar]
   Int_t           numnu;
   Int_t           mode;
   Int_t           ipnu[50];   //[numnu]
   Int_t           fqntwnd;
   Int_t           fqtwnd_iclstr[10];   //[fqntwnd]
   Int_t           fqtwnd_npeak[10];   //[fqntwnd]
   Float_t         fqtwnd_prftt0[10];   //[fqntwnd]
   Float_t         fqrcpar;
   Float_t         fqpi0par;
   Float_t         fqtwnd_prftpos[10][3];   //[fqntwnd]
   Float_t         fqtwnd[10][2];   //[fqntwnd]
   Float_t         fqtwnd_peakt0[10][10];   //[fqntwnd]
   Float_t         fqtwnd_peakiness[10][10];   //[fqntwnd]
   Int_t           fqnse;
   Int_t           fqitwnd[10];   //[fqnse]
   Int_t           fqipeak[10];   //[fqnse]
   Int_t           fqnhitpmt[10];   //[fqnse]
   Float_t         fqtotq[10];   //[fqnse]
   Float_t         fq0rtotmu[10];   //[fqnse]
   Float_t         fq0rnll[10];   //[fqnse]
   Int_t           fqn50[10];   //[fqnse]
   Float_t         fqq50[10];   //[fqnse]
   Int_t           fq1rpcflg[10][7];   //[fqnse]
   Float_t         fq1rmom[10][7];   //[fqnse]
   Float_t         fq1rt0[10][7];   //[fqnse]
   Float_t         fq1rtotmu[10][7];   //[fqnse]
   Float_t         fq1rnll[10][7];   //[fqnse]
   Float_t         fq1rpos[10][7][3];   //[fqnse]
   Float_t         fq1rdir[10][7][3];   //[fqnse]
   Float_t         fq1rdconv[10][7];   //[fqnse]
   Float_t         fq1reloss[10][7];   //[fqnse]
   Int_t           fqpi0pcflg[2];
   Float_t         fqpi0mom1[2];
   Float_t         fqpi0mom2[2];
   Float_t         fqpi0momtot[2];
   Float_t         fqpi0dconv1[2];
   Float_t         fqpi0dconv2[2];
   Float_t         fqpi0t0[2];
   Float_t         fqpi0totmu[2];
   Float_t         fqpi0nll[2];
   Float_t         fqpi0mass[2];
   Float_t         fqpi0photangle[2];
   Float_t         fqpi0pos[2][3];
   Float_t         fqpi0dir1[2][3];
   Float_t         fqpi0dir2[2][3];
   Float_t         fqpi0dirtot[2][3];
   Int_t           fqnmrfit;
   Int_t           fqmrifit[32];   //[fqnmrfit]
   Int_t           fqmrnring[32];   //[fqnmrfit]
   Int_t           fqmrpcflg[32];   //[fqnmrfit]
   Float_t         fqmrnll[32];   //[fqnmrfit]
   Float_t         fqmrtotmu[32];   //[fqnmrfit]
   Int_t           fqmrpid[32][6];   //[fqnmrfit]
   Float_t         fqmrmom[32][6];   //[fqnmrfit]
   Float_t         fqmrdconv[32][6];   //[fqnmrfit]
   Float_t         fqmreloss[32][6];   //[fqnmrfit]
   Float_t         fqmrt0[32][6];   //[fqnmrfit]
   Float_t         fqmrpos[32][6][3];   //[fqnmrfit]
   Float_t         fqmrdir[32][6][3];   //[fqnmrfit]
   Int_t           fqmsnfit;
   Int_t           fqmspcflg[5];   //[fqmsnfit]
   Int_t           fqmsnseg[5];   //[fqmsnfit]
   Int_t           fqmspid[5];   //[fqmsnfit]
   Int_t           fqmsifit[5];   //[fqmsnfit]
   Int_t           fqmsimer[5];   //[fqmsnfit]
   Float_t         fqmstotmu[5];   //[fqmsnfit]
   Float_t         fqmsnll[5];   //[fqmsnfit]
   Float_t         fqmsmom[5][20];   //[fqmsnfit]
   Float_t         fqmseloss[5][20];   //[fqmsnfit]
   Float_t         fqmst0[5][20];   //[fqmsnfit]
   Float_t         fqmspos[5][20][3];   //[fqmsnfit]
   Float_t         fqmsdir[5][20][3];   //[fqmsnfit]
   Float_t         attribute[1000];
   Int_t           ncomponent;
   Int_t           nsample;
   Int_t           nbin;
   Int_t           nvis;
   Int_t           nvmu;
   Int_t           nve;
   Int_t           nvgam;
   Int_t           nvpip;
   Int_t           nvpi0;
   Int_t           nvp;
   Int_t           nvk;
   Int_t           nvoth;
   Int_t           vispid[100];
   Float_t         fqwall;
   Float_t         fqtowall;
   Float_t         fq1rwall[10][7];
   Float_t         fq1rtowall[10][7];
   Float_t         fq1renu[2];
   Float_t         towallv[50];
   Float_t         wallv2;
   Float_t         evtweight;
   Int_t           best2RID;
   Int_t           passecut;
   Int_t           passe1rpicut;
   Int_t           passmucut;
   Float_t         fq1rperim[10][7];
   Float_t         fq1rmincone[10][7];
   // from Xiaoyue's skimmed tree with atmospheric weights
   Double_t        wgtosc1[4];
   Double_t        wgtflx[4];
   Float_t         oscpower[4];

   // List of branches
   TBranch        *b_nev;   //!
   TBranch        *b_nhitac;   //!
   TBranch        *b_npar;   //!
   TBranch        *b_wallv;   //!
   TBranch        *b_ipv;   //!
   TBranch        *b_posv;   //!
   TBranch        *b_dirv;   //!
   TBranch        *b_pmomv;   //!
   TBranch        *b_numnu;   //!
   TBranch        *b_mode;   //!
   TBranch        *b_ipnu;   //!
   TBranch        *b_fqntwnd;   //!
   TBranch        *b_fqtwnd_iclstr;   //!
   TBranch        *b_fqtwnd_npeak;   //!
   TBranch        *b_fqtwnd_prftt0;   //!
   TBranch        *b_fqtwnd_prftpos;   //!
   TBranch        *b_fqtwnd;   //!
   TBranch        *b_fqtwnd_peakt0;   //!
   TBranch        *b_fqtwnd_peakiness;   //!
   TBranch        *b_fqnse;   //!
   TBranch        *b_fqitwnd;   //!
   TBranch        *b_fqipeak;   //!
   TBranch        *b_fqnhitpmt;   //!
   TBranch        *b_fqtotq;   //!
   TBranch        *b_fq0rtotmu;   //!
   TBranch        *b_fq0rnll;   //!
   TBranch        *b_fqn50;   //!
   TBranch        *b_fqq50;   //!
   TBranch        *b_fq1rpcflg;   //!
   TBranch        *b_fq1rmom;   //!
   TBranch        *b_fq1rt0;   //!
   TBranch        *b_fq1rtotmu;   //!
   TBranch        *b_fq1rnll;   //!
   TBranch        *b_fq1rpos;   //!
   TBranch        *b_fq1rdir;   //!
   TBranch        *b_fq1rdconv;   //!
   TBranch        *b_fq1reloss;   //!
   TBranch        *b_fqpi0pcflg;   //!
   TBranch        *b_fqpi0mom1;   //!
   TBranch        *b_fqpi0mom2;   //!
   TBranch        *b_fqpi0momtot;   //!
   TBranch        *b_fqpi0dconv1;   //!
   TBranch        *b_fqpi0dconv2;   //!
   TBranch        *b_fqpi0t0;   //!
   TBranch        *b_fqpi0totmu;   //!
   TBranch        *b_fqpi0nll;   //!
   TBranch        *b_fqpi0mass;   //!
   TBranch        *b_fqpi0photangle;   //!
   TBranch        *b_fqpi0pos;   //!
   TBranch        *b_fqpi0dir1;   //!
   TBranch        *b_fqpi0dir2;   //!
   TBranch        *b_fqpi0dirtot;   //!
   TBranch        *b_fqnmrfit;   //!
   TBranch        *b_fqmrifit;   //!
   TBranch        *b_fqmrnring;   //!
   TBranch        *b_fqmrpcflg;   //!
   TBranch        *b_fqmrnll;   //!
   TBranch        *b_fqmrtotmu;   //!
   TBranch        *b_fqmrpid;   //!
   TBranch        *b_fqmrmom;   //!
   TBranch        *b_fqmrdconv;   //!
   TBranch        *b_fqmreloss;   //!
   TBranch        *b_fqmrt0;   //!
   TBranch        *b_fqmrpos;   //!
   TBranch        *b_fqmrdir;   //!
   TBranch        *b_fqmsnfit;   //!
   TBranch        *b_fqmspcflg;   //!
   TBranch        *b_fqmsnseg;   //!
   TBranch        *b_fqmspid;   //!
   TBranch        *b_fqmsifit;   //!
   TBranch        *b_fqmsimer;   //!
   TBranch        *b_fqmstotmu;   //!
   TBranch        *b_fqmsnll;   //!
   TBranch        *b_fqmsmom;   //!
   TBranch        *b_fqmseloss;   //!
   TBranch        *b_fqmst0;   //!
   TBranch        *b_fqmspos;   //!
   TBranch        *b_fqmsdir;   //!
   TBranch        *b_attribute;   //!
   TBranch        *b_ncomponent;   //!
   TBranch        *b_nsample;   //!
   TBranch        *b_nbin;   //!
   TBranch        *b_nvis;   //!
   TBranch        *b_nvmu;   //!
   TBranch        *b_nve;   //!
   TBranch        *b_nvgam;   //!
   TBranch        *b_nvpip;   //!
   TBranch        *b_nvpi0;   //!
   TBranch        *b_nvp;   //!
   TBranch        *b_nvk;   //!
   TBranch        *b_nvoth;   //!
   TBranch        *b_vispid;   //!
   TBranch        *b_fqwall;   //!
   TBranch        *b_fqtowall;   //!
   TBranch        *b_fq1rwall;   //!
   TBranch        *b_fq1rtowall;   //!
   TBranch        *b_towallv;   //!
   TBranch        *b_wallv2;   //!
   TBranch        *b_evtweight;   //!
   TBranch        *b_best2RID;   //!
   TBranch        *b_fq1rperim;   //!
   TBranch        *b_fq1rmincone;   //!

   fqProcessedEvent(TTree *tree=0);

   float getOscPower();
   int passMuCuts();
   int passECuts();
   void useImportantOnly();
   virtual ~fqProcessedEvent();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#ifdef CINTMODE
#include "fqProcessedEvent.cxx"
#endif


#endif // #ifdef fqProcessedEvent_cxx
