#define MRLoop_cxx
#include "MRLoop.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void MRLoop::fillHisto(int nbins, double xmin, double xmax, double ymin, double ymax){

  // make histos 
  hrc[0] = new TH2D("h1r","h1r",nbins,xmin,xmax,nbins,ymin,ymax);
  hrc[1] = new TH2D("hmr","hmr",nbins,xmin,xmax,nbins,ymin,ymax);
  hrcraw[0] = new TH2D("h1rraw","h1rraw",nbins,xmin,xmax,nbins,ymin,ymax);
  hrcraw[1] = new TH2D("hmrraw","hmrraw",nbins,xmin,xmax,nbins,ymin,ymax);

  // fill histos
  for (int i=0; i<fChain->GetEntries(); i++){
    fChain->GetEntry(i);
    if (nhitac>16) continue;
    if (fq1rmom[0][1]<30.) continue;
     double lnL = fq1rnll[0][2]-fqmrnll[best2RID];
     double pmin = TMath::Min(fqmrmom[best2RID][0],fqmrmom[best2RID][1]);
    if (nvis==1){
      hrc[0]->Fill(pmin,lnL);
      hrcraw[0]->Fill(pmin,lnL,evtweight);
    }
    if (nvis>1){
      hrc[1]->Fill(pmin,lnL,evtweight);
      hrcraw[1]->Fill(pmin,lnL);
    }
  }

  hrc[0]->Draw("colz");

  return;
}

void MRLoop::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L MRLoop.C
//      Root > MRLoop t
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
