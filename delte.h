//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Oct 14 11:30:08 2016 by ROOT version 5.28/00c
// from TTree MCMCpath/MCMCpath
// found on file: mcmctree.root
//////////////////////////////////////////////////////////

#ifndef delte_h
#define delte_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class delte {
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
   Int_t           parcomp[200];
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

   delte(TTree *tree=0);
   virtual ~delte();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef delte_cxx
delte::delte(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("mcmctree.root");
      if (!f) {
         f = new TFile("mcmctree.root");
      }
      tree = (TTree*)gDirectory->Get("MCMCpath");

   }
   Init(tree);
}

delte::~delte()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t delte::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t delte::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void delte::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("npars", &npars, &b_npars);
   fChain->SetBranchAddress("step", &step, &b_step);
   fChain->SetBranchAddress("par", par, &b_par);
   fChain->SetBranchAddress("pardiff", pardiff, &b_pardiff);
   fChain->SetBranchAddress("parnominal", parnominal, &b_parnominal);
   fChain->SetBranchAddress("parbin", parbin, &b_parbin);
   fChain->SetBranchAddress("parcomp", parcomp, &b_parcomp);
   fChain->SetBranchAddress("paratt", paratt, &b_paratt);
   fChain->SetBranchAddress("parindex", parindex, &b_parindex);
   fChain->SetBranchAddress("parsyst", parsyst, &b_parsyst);
   fChain->SetBranchAddress("logL", &logL, &b_logL);
   Notify();
}

Bool_t delte::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void delte::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t delte::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef delte_cxx
