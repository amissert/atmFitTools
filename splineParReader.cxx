#ifndef splineParReader_cxx
#define splineParReader_cxx


#include "splineParReader.h"


splineParReader::splineParReader(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("nom3_splineFactoryOut.root");
      if (!f) {
         f = new TFile("nom3_splineFactoryOut.root");
      }
      tree = (TTree*)gDirectory->Get("splinePars");

   }
   Init(tree);
}

splineParReader::~splineParReader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t splineParReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t splineParReader::LoadTree(Long64_t entry)
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

void splineParReader::Init(TTree *tree)
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

   fChain->SetBranchAddress("nbin", &nbin, &b_nbin);
   fChain->SetBranchAddress("nhistobins", &nhistobins, &b_nhistobins);
   fChain->SetBranchAddress("ncomponent", &ncomponent, &b_ncomponent);
   fChain->SetBranchAddress("nattribute", &nattribute, &b_nattribute);
   fChain->SetBranchAddress("nsample", &nsample, &b_nsample);
   if (fChain->GetListOfBranches()->FindObject("nmode"))
     fChain->SetBranchAddress("nmode", &nmode, &b_nmode);
   fChain->SetBranchAddress("nsystpar", &nsystpar, &b_nsystpar);
   fChain->SetBranchAddress("nsyspartot", &nsyspartot, &b_nsyspartot);
   fChain->SetBranchAddress("npoints", &npoints, &b_npoints);
   fChain->SetBranchAddress("systParValues", systParValues, &b_systParValues);
   fChain->SetBranchAddress("binWeight", binWeight, &b_binWeight);
   Notify();
}

Bool_t splineParReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void splineParReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

Int_t splineParReader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   entry++;
   return 1;
}



///////////////////////////////////////////////
// Draw 1D plot of the spline for a given bin
TH1D* splineParReader::drawSpline(int ibin){
  TH1D* h  = new TH1D();
  return h;
}


/////////////////////////////////////////////
// Draw a 2D plot of the histogram variation
TH2D* splineParReader::draw2D(int ievent){
 
  if (ievent>0) fChain->GetEntry(ievent);

  // setup histogram
  int nbinsx = nhistobins;
  double xmin   = 0;
  double xmax   = nhistobins;
  int nbinsy = npoints;
  double ymin = systParValues[0];
  double ymax = systParValues[npoints-1];
  TH2D* h2 = new TH2D("h2d","h2d",nbinsx,xmin,xmax,nbinsy,ymin,ymax);

  // fill histogram
  for (int iy=0; iy<npoints; iy++){
    for (int ix=1; ix<=nbinsx; ix++){
       double value = binWeight[iy][ix];
       h2->SetBinContent(ix,iy+1,value);
    }
  }

  // draw histogram
  h2->Draw("lego2");

  return h2;
}

void splineParReader::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L splineParReader.cxx
//      Root > splineParReader t
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
