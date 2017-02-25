#ifndef FQPROCEVENT_C
#define FQPROCEVENT_C

#include "fqProcessedEvent.h"

using namespace std;

/////////////////////////////////////////////////////
// Use only branches essential for cuts
void fqProcessedEvent::useImportantOnly(){
  
  fChain->SetBranchStatus("*",0);
  
  // fiTQun
  fChain->SetBranchStatus("fq1r*",1);
  fChain->SetBranchStatus("fqmrnring",1);

  // FC 
  fChain->SetBranchStatus("nhitac",1);

  // processed vars
  fChain->SetBranchStatus("n*",1);
  fChain->SetBranchStatus("at*",1);
  fChain->SetBranchStatus("ev*",1);
 
  return;
}

/////////////////////////////////////////////////////
// see if event passes numu cuts
int fqProcessedEvent::passMuCuts(){

   if (nhitac>16) return 0; //< FC cut
//   cout<<"pass 1"<<endl;
   if ((attribute[momentumIndex]) < 30.) return 0; //< Evis cut
//   cout<<"pass 2"<<endl;
   if (attribute[PIDIndex] > 0.) return 0; //< e/mu pid cut
//   cout<<"pass 3"<<endl;
   if (fqmrnring[0]>1) return 0; //< RC cut
//   cout<<"pass 3"<<endl;

   return 1;
}


/////////////////////////////////////////////////////
// see if event passes nue cuts
int fqProcessedEvent::passECuts(){

   if (nhitac>16) return 0; //< FC cut
   if ((attribute[momentumIndex]) < 100.) return 0; //< Evis cut
   if (attribute[PIDIndex] < 0.) return 0; //< e/mu pid cut
   if (fqmrnring[0]>1) return 0; //< RC cut
   if (fqpi0par>0) return 0;
//   float Lpi0 = fq1rnll[0][1] - fqpi0nll[0]; 
//   float pmome = fq1rmom[0][1];
//   float pi0par = Lpi0 - 70. - ((140.-70.)/(40.- 120.))*(pmome - 120.);

   return 1;
}



///////////////////////////////////////////////
// Get derivative weight (working on it)
float fqProcessedEvent::getOscPower(){

  if (TMath::Abs(mode)==1) return 1.;
  return 0.;

}

fqProcessedEvent::fqProcessedEvent(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("t2kmc_numu_ppmc_2_.root");
      if (!f) {
         f = new TFile("t2kmc_numu_ppmc_2_.root");
      }
      tree = (TTree*)gDirectory->Get("h1");

   }
   Init(tree);
}

fqProcessedEvent::~fqProcessedEvent()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t fqProcessedEvent::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t fqProcessedEvent::LoadTree(Long64_t entry)
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

void fqProcessedEvent::Init(TTree *tree)
{


   momentumIndex=3;
   PIDIndex=0;

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
   fChain->SetBranchAddress("fqrcpar",&fqrcpar);
   fChain->SetBranchAddress("fqpi0par",&fqpi0par);
   fChain->SetBranchAddress("nev", &nev, &b_nev);
   fChain->SetBranchAddress("nhitac", &nhitac, &b_nhitac);
   fChain->SetBranchAddress("npar", &npar, &b_npar);
   fChain->SetBranchAddress("wallv", &wallv, &b_wallv);
   fChain->SetBranchAddress("ipv", ipv, &b_ipv);
   fChain->SetBranchAddress("posv", posv, &b_posv);
   fChain->SetBranchAddress("dirv", dirv, &b_dirv);
   fChain->SetBranchAddress("pmomv", pmomv, &b_pmomv);
   fChain->SetBranchAddress("numnu", &numnu, &b_numnu);
   fChain->SetBranchAddress("mode", &mode, &b_mode);
   fChain->SetBranchAddress("ipnu", ipnu, &b_ipnu);
   fChain->SetBranchAddress("fqntwnd", &fqntwnd, &b_fqntwnd);
   fChain->SetBranchAddress("fqtwnd_iclstr", fqtwnd_iclstr, &b_fqtwnd_iclstr);
   fChain->SetBranchAddress("fqtwnd_npeak", fqtwnd_npeak, &b_fqtwnd_npeak);
   fChain->SetBranchAddress("fqtwnd_prftt0", fqtwnd_prftt0, &b_fqtwnd_prftt0);
   fChain->SetBranchAddress("fqtwnd_prftpos", fqtwnd_prftpos, &b_fqtwnd_prftpos);
   fChain->SetBranchAddress("fqtwnd", fqtwnd, &b_fqtwnd);
   fChain->SetBranchAddress("fqtwnd_peakt0", fqtwnd_peakt0, &b_fqtwnd_peakt0);
   fChain->SetBranchAddress("fqtwnd_peakiness", fqtwnd_peakiness, &b_fqtwnd_peakiness);
   fChain->SetBranchAddress("fqnse", &fqnse, &b_fqnse);
   fChain->SetBranchAddress("fqitwnd", fqitwnd, &b_fqitwnd);
   fChain->SetBranchAddress("fqipeak", fqipeak, &b_fqipeak);
   fChain->SetBranchAddress("fqnhitpmt", fqnhitpmt, &b_fqnhitpmt);
   fChain->SetBranchAddress("fqtotq", fqtotq, &b_fqtotq);
   fChain->SetBranchAddress("fq0rtotmu", fq0rtotmu, &b_fq0rtotmu);
   fChain->SetBranchAddress("fq0rnll", fq0rnll, &b_fq0rnll);
   fChain->SetBranchAddress("fqn50", fqn50, &b_fqn50);
   fChain->SetBranchAddress("fqq50", fqq50, &b_fqq50);
   fChain->SetBranchAddress("fq1rpcflg", fq1rpcflg, &b_fq1rpcflg);
   fChain->SetBranchAddress("fq1rmom", fq1rmom, &b_fq1rmom);
   fChain->SetBranchAddress("fq1rt0", fq1rt0, &b_fq1rt0);
   fChain->SetBranchAddress("fq1rtotmu", fq1rtotmu, &b_fq1rtotmu);
   fChain->SetBranchAddress("fq1rnll", fq1rnll, &b_fq1rnll);
   fChain->SetBranchAddress("fq1rpos", fq1rpos, &b_fq1rpos);
   fChain->SetBranchAddress("fq1renu", fq1renu);
   fChain->SetBranchAddress("passecut", &passecut);
   fChain->SetBranchAddress("passe1rpicut", &passe1rpicut);
   fChain->SetBranchAddress("passmucut", &passmucut);
   fChain->SetBranchAddress("fq1rdir", fq1rdir, &b_fq1rdir);
   fChain->SetBranchAddress("fq1rdconv", fq1rdconv, &b_fq1rdconv);
   fChain->SetBranchAddress("fq1reloss", fq1reloss, &b_fq1reloss);
   fChain->SetBranchAddress("fqpi0pcflg", fqpi0pcflg, &b_fqpi0pcflg);
   fChain->SetBranchAddress("fqpi0mom1", fqpi0mom1, &b_fqpi0mom1);
   fChain->SetBranchAddress("fqpi0mom2", fqpi0mom2, &b_fqpi0mom2);
   fChain->SetBranchAddress("fqpi0momtot", fqpi0momtot, &b_fqpi0momtot);
   fChain->SetBranchAddress("fqpi0dconv1", fqpi0dconv1, &b_fqpi0dconv1);
   fChain->SetBranchAddress("fqpi0dconv2", fqpi0dconv2, &b_fqpi0dconv2);
   fChain->SetBranchAddress("fqpi0t0", fqpi0t0, &b_fqpi0t0);
   fChain->SetBranchAddress("fqpi0totmu", fqpi0totmu, &b_fqpi0totmu);
   fChain->SetBranchAddress("fqpi0nll", fqpi0nll, &b_fqpi0nll);
   fChain->SetBranchAddress("fqpi0mass", fqpi0mass, &b_fqpi0mass);
   fChain->SetBranchAddress("fqpi0photangle", fqpi0photangle, &b_fqpi0photangle);
   fChain->SetBranchAddress("fqpi0pos", fqpi0pos, &b_fqpi0pos);
   fChain->SetBranchAddress("fqpi0dir1", fqpi0dir1, &b_fqpi0dir1);
   fChain->SetBranchAddress("fqpi0dir2", fqpi0dir2, &b_fqpi0dir2);
   fChain->SetBranchAddress("fqpi0dirtot", fqpi0dirtot, &b_fqpi0dirtot);
   fChain->SetBranchAddress("fqnmrfit", &fqnmrfit, &b_fqnmrfit);
//   fChain->SetBranchAddress("fqmrifit", fqmrifit, &b_fqmrifit);
   fChain->SetBranchAddress("fqmrnring", fqmrnring, &b_fqmrnring);
   fChain->SetBranchAddress("fqmrpcflg", fqmrpcflg, &b_fqmrpcflg);
   fChain->SetBranchAddress("fqmrnll", fqmrnll, &b_fqmrnll);
   fChain->SetBranchAddress("fqmrtotmu", fqmrtotmu, &b_fqmrtotmu);
   fChain->SetBranchAddress("fqmrpid", fqmrpid, &b_fqmrpid);
   fChain->SetBranchAddress("fqmrmom", fqmrmom, &b_fqmrmom);
   fChain->SetBranchAddress("fqmrdconv", fqmrdconv, &b_fqmrdconv);
   fChain->SetBranchAddress("fqmreloss", fqmreloss, &b_fqmreloss);
   fChain->SetBranchAddress("fqmrt0", fqmrt0, &b_fqmrt0);
   fChain->SetBranchAddress("fqmrpos", fqmrpos, &b_fqmrpos);
   fChain->SetBranchAddress("fqmrdir", fqmrdir, &b_fqmrdir);
   fChain->SetBranchAddress("fqmsnfit", &fqmsnfit, &b_fqmsnfit);
   fChain->SetBranchAddress("fqmspcflg", fqmspcflg, &b_fqmspcflg);
   fChain->SetBranchAddress("fqmsnseg", fqmsnseg, &b_fqmsnseg);
   fChain->SetBranchAddress("fqmspid", fqmspid, &b_fqmspid);
   fChain->SetBranchAddress("fqmsifit", fqmsifit, &b_fqmsifit);
   fChain->SetBranchAddress("fqmsimer", fqmsimer, &b_fqmsimer);
   fChain->SetBranchAddress("fqmstotmu", fqmstotmu, &b_fqmstotmu);
   fChain->SetBranchAddress("fqmsnll", fqmsnll, &b_fqmsnll);
   fChain->SetBranchAddress("fqmsmom", fqmsmom, &b_fqmsmom);
   fChain->SetBranchAddress("fqmseloss", fqmseloss, &b_fqmseloss);
   fChain->SetBranchAddress("fqmst0", fqmst0, &b_fqmst0);
   fChain->SetBranchAddress("fqmspos", fqmspos, &b_fqmspos);
   fChain->SetBranchAddress("fqmsdir", fqmsdir, &b_fqmsdir);
   fChain->SetBranchAddress("attribute", attribute, &b_attribute);
   fChain->SetBranchAddress("ncomponent", &ncomponent, &b_ncomponent);
   fChain->SetBranchAddress("nsample", &nsample, &b_nsample);
   fChain->SetBranchAddress("nbin", &nbin, &b_nbin);
   fChain->SetBranchAddress("nvis", &nvis, &b_nvis);
   fChain->SetBranchAddress("nvmu", &nvmu, &b_nvmu);
   fChain->SetBranchAddress("nve", &nve, &b_nve);
   fChain->SetBranchAddress("nvgam", &nvgam, &b_nvgam);
   fChain->SetBranchAddress("nvpip", &nvpip, &b_nvpip);
   fChain->SetBranchAddress("nvpi0", &nvpi0, &b_nvpi0);
   fChain->SetBranchAddress("nvp", &nvp, &b_nvp);
   fChain->SetBranchAddress("nvk", &nvk, &b_nvk);
   //fChain->SetBranchAddress("nvoth", &nvoth, &b_nvoth);;
//   fChain->SetBranchAddress("vispid", vispid, &b_vispid);
   fChain->SetBranchAddress("fqwall", &fqwall, &b_fqwall);
   fChain->SetBranchAddress("fqtowall", &fqtowall, &b_fqtowall);
   fChain->SetBranchAddress("fq1rwall", fq1rwall, &b_fq1rwall);
   fChain->SetBranchAddress("fq1rtowall", fq1rtowall, &b_fq1rtowall);
   fChain->SetBranchAddress("towallv", towallv, &b_towallv);
   fChain->SetBranchAddress("wallv2", &wallv2, &b_wallv2);
   fChain->SetBranchAddress("evtweight", &evtweight, &b_evtweight);
   fChain->SetBranchAddress("best2RID", &best2RID, &b_best2RID);
   fChain->SetBranchAddress("oscpower",oscpower);
//   fChain->SetBranchAddress("fq1rperim", fq1rperim, &b_fq1rperim);
//   fChain->SetBranchAddress("fq1rmincone", fq1rmincone, &b_fq1rmincone);
//#ifdef USE_ATM_WEIGHTS
//   fChain->SetBranchAddress("wgtosc1", wgtosc1);
//   fChain->SetBranchAddress("wgtflx", wgtflx);
//#endif
   Notify();
}

Bool_t fqProcessedEvent::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void fqProcessedEvent::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

Int_t fqProcessedEvent::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   entry++;
   return 1;
}



#endif
