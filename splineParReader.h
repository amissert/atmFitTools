//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Dec  8 18:55:13 2015 by ROOT version 5.28/00c
// from TTree splinePars/spinePars
// found on file: nom3_splineFactoryOut.root
//////////////////////////////////////////////////////////

#ifndef splineParReader_h
#define splineParReader_h
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "shared.h"
#include "TH2D.h"
#include "TH1D.h"

class splineParReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           nbin;
   Int_t           nhistobins;
   Int_t           ncomponent;
   Int_t           nattribute;
   Int_t           nsample;
   Int_t           nmode;
   Int_t           nsystpar;
   Int_t           nsyspartot;
   Int_t           npoints;
//#ifndef T2K
//   Double_t        systParValues[5];
//   Double_t        binWeight[5][300];
//#else
   Double_t        systParValues[NPTSMAX];
   Double_t        binWeight[NPTSMAX][NHBINSMAX];
//#endif


   // List of branches
   TBranch        *b_nbin;   //!
   TBranch        *b_nhistobins;   //!
   TBranch        *b_ncomponent;   //!
   TBranch        *b_nattribute;   //!
   TBranch        *b_nsample;   //!
   TBranch        *b_nmode;
   TBranch        *b_nsystpar;   //!
   TBranch        *b_nsyspartot;   //!
   TBranch        *b_npoints;   //!
   TBranch        *b_systParValues;   //!
   TBranch        *b_binWeight;   //!

   splineParReader(TTree *tree=0);
   virtual ~splineParReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     Loop();


   // plotting functions for debugging etc.
   TH2D* draw2D(int ientry = -1);
   TH1D* drawSpline(int ibin);
};

#ifdef CINTMODE
#include "splineParReader.cxx"
#endif

#endif //
