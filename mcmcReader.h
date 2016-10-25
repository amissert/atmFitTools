//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct 19 14:22:33 2016 by ROOT version 5.28/00c
// from TTree MCMCpath/MCMCpath
// found on file: fullfit_allpars_norm_nosyst.root
//////////////////////////////////////////////////////////

#ifndef mcmcReader_h
#define mcmcReader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

class mcmcReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           npars;
   Int_t           step;
   Float_t         par[500];
   Float_t         pardiff[500];
   Float_t         parnominal[500];
   Int_t           parbin[500];
   Int_t           parcomp[500];
   Int_t           paratt[500];
   Int_t           parindex[500];
   Int_t           parsyst[500];
   Float_t         logL;

   // List of branches
   TBranch        *b_npars;   //!
   TBranch        *b_step;   //!
   TBranch        *b_par;   //!
   TBranch        *b_pardiff;   //!
   TBranch        *b_parnominal;   //!
   TBranch        *b_parbin;   //!
   TBranch        *b_parcomp;   //!
   TBranch        *b_paratt;   //!
   TBranch        *b_parindex;   //!
   TBranch        *b_parsyst;   //!
   TBranch        *b_logL;   //!

   mcmcReader(TTree *tree=0);
   virtual ~mcmcReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#ifdef CINTMODE
#include "mcmcReader.cxx"
#endif


#endif // #ifdef mcmcReader_cxx
