#ifndef mcmcReader_cxx
#define mcmcReader_cxx

#include "mcmcReader.h"


mcmcReader::mcmcReader(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("fullfit_allpars_norm_nosyst.root");
      if (!f) {
         f = new TFile("fullfit_allpars_norm_nosyst.root");
      }
      tree = (TTree*)gDirectory->Get("MCMCpath");

   }
   Init(tree);
}

mcmcReader::~mcmcReader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t mcmcReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t mcmcReader::LoadTree(Long64_t entry)
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

void mcmcReader::Init(TTree *tree)
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

Bool_t mcmcReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void mcmcReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t mcmcReader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
;

void mcmcReader::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L mcmcReader.C
//      Root > mcmcReader t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}



#endif
