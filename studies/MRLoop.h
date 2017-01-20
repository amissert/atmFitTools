//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Nov 13 16:31:42 2016 by ROOT version 5.28/00c
// from TChain h1/
//////////////////////////////////////////////////////////

#ifndef MRLoop_h
#define MRLoop_h

#include <TROOT.h>
#include <TChain.h>
#include <TH2D.h>
#include <TFile.h>

class MRLoop {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           nring;
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
   Int_t           nscndprt;
   Int_t           itrkscnd[1000];   //[nscndprt]
   Int_t           istakscnd[1000];   //[nscndprt]
   Float_t         vtxscnd[1000][3];   //[nscndprt]
   Float_t         pscnd[1000][3];   //[nscndprt]
   Int_t           iprtscnd[1000];   //[nscndprt]
   Float_t         tscnd[1000];   //[nscndprt]
   Int_t           iprntprt[1000];   //[nscndprt]
   Int_t           lmecscnd[1000];   //[nscndprt]
   Int_t           iprnttrk[1000];   //[nscndprt]
   Int_t           iorgprt[1000];   //[nscndprt]
   Int_t           iflgscnd[1000];   //[nscndprt]
   Int_t           iprntidx[1000];   //[nscndprt]
   Int_t           nchilds[1000];   //[nscndprt]
   Int_t           ichildidx[1000];   //[nscndprt]
   Int_t           Npvcscnd;
   Int_t           Nchildsvc[100];   //[Npvcscnd]
   Int_t           Ichildidxvc[100];   //[Npvcscnd]
   Int_t           Npvc;
   Int_t           Ipvc[100];   //[Npvc]
   Int_t           Ichvc[100];   //[Npvc]
   Int_t           Iorgvc[100];   //[Npvc]
   Int_t           Iflvc[100];   //[Npvc]
   Float_t         Abspvc[100];   //[Npvc]
   Float_t         Pvc[100][3];   //[Npvc]
   Int_t           fqntwnd;
   Int_t           fqtwnd_iclstr[10];   //[fqntwnd]
   Int_t           fqtwnd_npeak[10];   //[fqntwnd]
   Float_t         fqtwnd_prftt0[10];   //[fqntwnd]
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
   Double_t        totwgt;
   Double_t        oscwgt;
   Double_t        normwgt;
   Float_t         attribute[100];
   Float_t         fqrcpar;
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
   Double_t        visbrightness[17];   //[nvis]
   Double_t        viswall[100];
   Double_t        vistowall[17];   //[nvis]
   Double_t        vismrwall1;
   Double_t        vismrwall2;
   Double_t        vismrtowall1;
   Double_t        vismrtowall2;
   Double_t        vismrwallmin;
   Double_t        vismrtowallmin;
   Int_t           nvisscnd;
   Double_t        vismrbrightness;
   Int_t           vismrpid1;
   Int_t           vismrpid2;
   Double_t        vismrt1;
   Double_t        vismrt2;
   Int_t           vismrtype1;
   Int_t           vismrtype2;
   Double_t        vistime[17];   //[nvis]
   Int_t           vispid[17];   //[nvis]
   Int_t           visscndpid[17];   //[nvis]
   Int_t           visscndparentid[17];   //[nvis]
   Float_t         fqwall;
   Float_t         fqtowall;
   Float_t         fq1rwall[10][7];
   Float_t         fq1rtowall[10][7];
   Float_t         fq1renu[2];
   Float_t         towallv[50];
   Float_t         wallv2;
   Float_t         evtweight;
   Int_t           best2RID;
   Float_t         fq1rperim[10][7];
   Float_t         fq1rmincone[10][7];

   // List of branches
   TBranch        *b_nring;   //!
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
   TBranch        *b_nscndprt;   //!
   TBranch        *b_itrkscnd;   //!
   TBranch        *b_istakscnd;   //!
   TBranch        *b_vtxscnd;   //!
   TBranch        *b_pscnd;   //!
   TBranch        *b_iprtscnd;   //!
   TBranch        *b_tscnd;   //!
   TBranch        *b_iprntprt;   //!
   TBranch        *b_lmecscnd;   //!
   TBranch        *b_iprnttrk;   //!
   TBranch        *b_iorgprt;   //!
   TBranch        *b_iflgscnd;   //!
   TBranch        *b_iprntidx;   //!
   TBranch        *b_nchilds;   //!
   TBranch        *b_ichildidx;   //!
   TBranch        *b_Npvcscnd;   //!
   TBranch        *b_Nchildsvc;   //!
   TBranch        *b_Ichildidxvc;   //!
   TBranch        *b_Npvc;   //!
   TBranch        *b_Ipvc;   //!
   TBranch        *b_Ichvc;   //!
   TBranch        *b_Iorgvc;   //!
   TBranch        *b_Iflvc;   //!
   TBranch        *b_Abspvc;   //!
   TBranch        *b_Pvc;   //!
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
   TBranch        *b_totwgt;   //!
   TBranch        *b_oscwgt;   //!
   TBranch        *b_normwgt;   //!
   TBranch        *b_attribute;   //!
   TBranch        *b_fqrcpar;   //!
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
   TBranch        *b_visbrightness;   //!
   TBranch        *b_viswall;   //!
   TBranch        *b_vistowall;   //!
   TBranch        *b_vismrwall1;   //!
   TBranch        *b_vismrwall2;   //!
   TBranch        *b_vismrtowall1;   //!
   TBranch        *b_vismrtowall2;   //!
   TBranch        *b_vismrwallmin;   //!
   TBranch        *b_vismrtowallmin;   //!
   TBranch        *b_nvisscnd;   //!
   TBranch        *b_vismrbrightness;   //!
   TBranch        *b_vismrpid1;   //!
   TBranch        *b_vismrpid2;   //!
   TBranch        *b_vismrt1;   //!
   TBranch        *b_vismrt2;   //!
   TBranch        *b_vismrtype1;   //!
   TBranch        *b_vismrtype2;   //!
   TBranch        *b_vistime;   //!
   TBranch        *b_vispid;   //!
   TBranch        *b_visscndpid;   //!
   TBranch        *b_visscndparentid;   //!
   TBranch        *b_fqwall;   //!
   TBranch        *b_fqtowall;   //!
   TBranch        *b_fq1rwall;   //!
   TBranch        *b_fq1rtowall;   //!
   TBranch        *b_fq1renu;   //!
   TBranch        *b_towallv;   //!
   TBranch        *b_wallv2;   //!
   TBranch        *b_evtweight;   //!
   TBranch        *b_best2RID;   //!
   TBranch        *b_fq1rperim;   //!
   TBranch        *b_fq1rmincone;   //!

   TH2D* hrc[2]; //< weighted histograms
   TH2D* hrcraw[2]; //< statistical histograms 
   MRLoop(TTree *tree=0);
   virtual ~MRLoop();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   void fillHisto(int nbins, double xmin, double xmax, double ymin, double ymax);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MRLoop_cxx
MRLoop::MRLoop(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/nfs/data41/t2k/amissert/t2kmc/processed/nov10/t2k_reduced_ppmc_917.root");
      if (!f) {
         f = new TFile("/nfs/data41/t2k/amissert/t2kmc/processed/nov10/t2k_reduced_ppmc_917.root");
      }
      tree = (TTree*)gDirectory->Get("h1");

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("h1","");
      chain->Add("/nfs/data41/t2k/amissert/t2kmc/processed/nov10/t2k_reduced_ppmc_1017.root/h1");
      chain->Add("/nfs/data41/t2k/amissert/t2kmc/processed/nov10/t2k_reduced_ppmc_1117.root/h1");
      chain->Add("/nfs/data41/t2k/amissert/t2kmc/processed/nov10/t2k_reduced_ppmc_117.root/h1");
      chain->Add("/nfs/data41/t2k/amissert/t2kmc/processed/nov10/t2k_reduced_ppmc_1217.root/h1");
      chain->Add("/nfs/data41/t2k/amissert/t2kmc/processed/nov10/t2k_reduced_ppmc_1317.root/h1");
      chain->Add("/nfs/data41/t2k/amissert/t2kmc/processed/nov10/t2k_reduced_ppmc_1417.root/h1");
      chain->Add("/nfs/data41/t2k/amissert/t2kmc/processed/nov10/t2k_reduced_ppmc_1517.root/h1");
      chain->Add("/nfs/data41/t2k/amissert/t2kmc/processed/nov10/t2k_reduced_ppmc_1617.root/h1");
      chain->Add("/nfs/data41/t2k/amissert/t2kmc/processed/nov10/t2k_reduced_ppmc_17.root/h1");
      chain->Add("/nfs/data41/t2k/amissert/t2kmc/processed/nov10/t2k_reduced_ppmc_1717.root/h1");
      chain->Add("/nfs/data41/t2k/amissert/t2kmc/processed/nov10/t2k_reduced_ppmc_1817.root/h1");
      chain->Add("/nfs/data41/t2k/amissert/t2kmc/processed/nov10/t2k_reduced_ppmc_1917.root/h1");
      chain->Add("/nfs/data41/t2k/amissert/t2kmc/processed/nov10/t2k_reduced_ppmc_2017.root/h1");
      chain->Add("/nfs/data41/t2k/amissert/t2kmc/processed/nov10/t2k_reduced_ppmc_2117.root/h1");
      chain->Add("/nfs/data41/t2k/amissert/t2kmc/processed/nov10/t2k_reduced_ppmc_217.root/h1");
      chain->Add("/nfs/data41/t2k/amissert/t2kmc/processed/nov10/t2k_reduced_ppmc_2217.root/h1");
      chain->Add("/nfs/data41/t2k/amissert/t2kmc/processed/nov10/t2k_reduced_ppmc_2317.root/h1");
      chain->Add("/nfs/data41/t2k/amissert/t2kmc/processed/nov10/t2k_reduced_ppmc_2417.root/h1");
      chain->Add("/nfs/data41/t2k/amissert/t2kmc/processed/nov10/t2k_reduced_ppmc_2517.root/h1");
      chain->Add("/nfs/data41/t2k/amissert/t2kmc/processed/nov10/t2k_reduced_ppmc_2617.root/h1");
      chain->Add("/nfs/data41/t2k/amissert/t2kmc/processed/nov10/t2k_reduced_ppmc_2717.root/h1");
      chain->Add("/nfs/data41/t2k/amissert/t2kmc/processed/nov10/t2k_reduced_ppmc_317.root/h1");
      chain->Add("/nfs/data41/t2k/amissert/t2kmc/processed/nov10/t2k_reduced_ppmc_417.root/h1");
      chain->Add("/nfs/data41/t2k/amissert/t2kmc/processed/nov10/t2k_reduced_ppmc_517.root/h1");
      chain->Add("/nfs/data41/t2k/amissert/t2kmc/processed/nov10/t2k_reduced_ppmc_617.root/h1");
      chain->Add("/nfs/data41/t2k/amissert/t2kmc/processed/nov10/t2k_reduced_ppmc_717.root/h1");
      chain->Add("/nfs/data41/t2k/amissert/t2kmc/processed/nov10/t2k_reduced_ppmc_817.root/h1");
      chain->Add("/nfs/data41/t2k/amissert/t2kmc/processed/nov10/t2k_reduced_ppmc_917.root/h1");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

MRLoop::~MRLoop()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MRLoop::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MRLoop::LoadTree(Long64_t entry)
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

void MRLoop::Init(TTree *tree)
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

   fChain->SetBranchAddress("nring", &nring, &b_nring);
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
   fChain->SetBranchAddress("nscndprt", &nscndprt, &b_nscndprt);
   fChain->SetBranchAddress("itrkscnd", itrkscnd, &b_itrkscnd);
   fChain->SetBranchAddress("istakscnd", istakscnd, &b_istakscnd);
   fChain->SetBranchAddress("vtxscnd", vtxscnd, &b_vtxscnd);
   fChain->SetBranchAddress("pscnd", pscnd, &b_pscnd);
   fChain->SetBranchAddress("iprtscnd", iprtscnd, &b_iprtscnd);
   fChain->SetBranchAddress("tscnd", tscnd, &b_tscnd);
   fChain->SetBranchAddress("iprntprt", iprntprt, &b_iprntprt);
   fChain->SetBranchAddress("lmecscnd", lmecscnd, &b_lmecscnd);
   fChain->SetBranchAddress("iprnttrk", iprnttrk, &b_iprnttrk);
   fChain->SetBranchAddress("iorgprt", iorgprt, &b_iorgprt);
   fChain->SetBranchAddress("iflgscnd", iflgscnd, &b_iflgscnd);
   fChain->SetBranchAddress("iprntidx", iprntidx, &b_iprntidx);
   fChain->SetBranchAddress("nchilds", nchilds, &b_nchilds);
   fChain->SetBranchAddress("ichildidx", ichildidx, &b_ichildidx);
   fChain->SetBranchAddress("Npvcscnd", &Npvcscnd, &b_Npvcscnd);
   fChain->SetBranchAddress("Nchildsvc", Nchildsvc, &b_Nchildsvc);
   fChain->SetBranchAddress("Ichildidxvc", Ichildidxvc, &b_Ichildidxvc);
   fChain->SetBranchAddress("Npvc", &Npvc, &b_Npvc);
   fChain->SetBranchAddress("Ipvc", Ipvc, &b_Ipvc);
   fChain->SetBranchAddress("Ichvc", Ichvc, &b_Ichvc);
   fChain->SetBranchAddress("Iorgvc", Iorgvc, &b_Iorgvc);
   fChain->SetBranchAddress("Iflvc", Iflvc, &b_Iflvc);
   fChain->SetBranchAddress("Abspvc", Abspvc, &b_Abspvc);
   fChain->SetBranchAddress("Pvc", Pvc, &b_Pvc);
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
   fChain->SetBranchAddress("fqmrifit", fqmrifit, &b_fqmrifit);
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
   fChain->SetBranchAddress("totwgt", &totwgt, &b_totwgt);
   fChain->SetBranchAddress("oscwgt", &oscwgt, &b_oscwgt);
   fChain->SetBranchAddress("normwgt", &normwgt, &b_normwgt);
   fChain->SetBranchAddress("attribute", attribute, &b_attribute);
   fChain->SetBranchAddress("fqrcpar", &fqrcpar, &b_fqrcpar);
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
   fChain->SetBranchAddress("visbrightness", visbrightness, &b_visbrightness);
   fChain->SetBranchAddress("viswall", viswall, &b_viswall);
   fChain->SetBranchAddress("vistowall", vistowall, &b_vistowall);
   fChain->SetBranchAddress("vismrwall1", &vismrwall1, &b_vismrwall1);
   fChain->SetBranchAddress("vismrwall2", &vismrwall2, &b_vismrwall2);
   fChain->SetBranchAddress("vismrtowall1", &vismrtowall1, &b_vismrtowall1);
   fChain->SetBranchAddress("vismrtowall2", &vismrtowall2, &b_vismrtowall2);
   fChain->SetBranchAddress("vismrwallmin", &vismrwallmin, &b_vismrwallmin);
   fChain->SetBranchAddress("vismrtowallmin", &vismrtowallmin, &b_vismrtowallmin);
   fChain->SetBranchAddress("nvisscnd", &nvisscnd, &b_nvisscnd);
   fChain->SetBranchAddress("vismrbrightness", &vismrbrightness, &b_vismrbrightness);
   fChain->SetBranchAddress("vismrpid1", &vismrpid1, &b_vismrpid1);
   fChain->SetBranchAddress("vismrpid2", &vismrpid2, &b_vismrpid2);
   fChain->SetBranchAddress("vismrt1", &vismrt1, &b_vismrt1);
   fChain->SetBranchAddress("vismrt2", &vismrt2, &b_vismrt2);
   fChain->SetBranchAddress("vismrtype1", &vismrtype1, &b_vismrtype1);
   fChain->SetBranchAddress("vismrtype2", &vismrtype2, &b_vismrtype2);
   fChain->SetBranchAddress("vistime", vistime, &b_vistime);
   fChain->SetBranchAddress("vispid", vispid, &b_vispid);
   fChain->SetBranchAddress("visscndpid", visscndpid, &b_visscndpid);
   fChain->SetBranchAddress("visscndparentid", visscndparentid, &b_visscndparentid);
   fChain->SetBranchAddress("fqwall", &fqwall, &b_fqwall);
   fChain->SetBranchAddress("fqtowall", &fqtowall, &b_fqtowall);
   fChain->SetBranchAddress("fq1rwall", fq1rwall, &b_fq1rwall);
   fChain->SetBranchAddress("fq1rtowall", fq1rtowall, &b_fq1rtowall);
   fChain->SetBranchAddress("fq1renu", fq1renu, &b_fq1renu);
   fChain->SetBranchAddress("towallv", towallv, &b_towallv);
   fChain->SetBranchAddress("wallv2", &wallv2, &b_wallv2);
   fChain->SetBranchAddress("evtweight", &evtweight, &b_evtweight);
   fChain->SetBranchAddress("best2RID", &best2RID, &b_best2RID);
   fChain->SetBranchAddress("fq1rperim", fq1rperim, &b_fq1rperim);
   fChain->SetBranchAddress("fq1rmincone", fq1rmincone, &b_fq1rmincone);
   Notify();
}

Bool_t MRLoop::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MRLoop::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MRLoop::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MRLoop_cxx
