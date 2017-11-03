{

gROOT->ProcessLine(".L histoCompare.cxx+");
gROOT->ProcessLine(".L atmFitPars.cxx+");
gROOT->ProcessLine(".L histoManager.cxx+");
gROOT->ProcessLine(".L modHistoArrayFV.cxx+");
gROOT->ProcessLine(".L mcmcApply.cxx+");
gROOT->ProcessLine(".L TH2FV.cxx+");
gROOT->ProcessLine(".L SKError.cxx+");
gROOT->ProcessLine(".L toyMC.cxx+");
gROOT->ProcessLine(".x ~/style.c");
gROOT->ProcessLine(".x ~/corr.c");


////////////////////////////////////////////////////////////////////////////////////////////////////////
// Style
//gStyle->SetTitleW(0.85);
//gStyle->SetTitleH(0.08);
//gStyle->SetOptTitle(1);


////////////////////////////////////////////////////////////////////////////////////////////////////////
// parameters //
TString AtmMCFiles = "/nfs/data41/t2k/amissert/processed/wetrun/jun27newrc_allevis/*ppmc*.root";
//TString MCMCFiles = "./run/results/wetrun_lowevis/mcmc_summary.root";
TString MCMCFiles = "./run/results/wetrun_newrc/demcmc_summary_run3.root";
int   NMCEvents = 1e9;
int   NMCMCPoints = 1700;
int   index_of_pidpar = 0;
int   index_of_pi0par = 1;
int   index_of_pippar = 2;
int   index_of_rcpar  = 3;
int   flgRCWeight = 1;
int   bindef = 2; //< 1 for tn-186, 2 to include NC bins
////////////////////////////////////////////////////////////////////////////////////////////////////////


// run toyMC
toyMC* toy = new toyMC();
toy->indexPIDPar = index_of_pidpar;
toy->indexPi0Par = index_of_pi0par;
toy->indexPiPPar = index_of_pippar;
toy->indexRCPar  = index_of_rcpar;

// mc files
TChain* mcfiles = new TChain("h1");
mcfiles->Add(AtmMCFiles.Data());

// mcmc pars
TChain* parfiles = new TChain("MCMCpath");
parfiles->Add(MCMCFiles.Data());
toy->setChains(mcfiles,parfiles,NMCEvents);
toy->setAtmFitPars("wetrun_newrc.dat");
mcfiles->GetEntry(500); // initialize some parameters
parfiles->GetEntry(500); // initialize some parameters


//toy->modifier->setAllAlphaFlags(false);
//toy->modifier->setAttFlgs(0,false);
//toy->modifier->setAttFlgs(1,false);
//toy->modifier->setAttFlgs(2,false);
toy->flgRCWeight = flgRCWeight;

toy->fillSKErrors(NMCMCPoints,bindef,0,0);

}
