{
gROOT->ProcessLine(".L histoCompare.cxx+");
gROOT->ProcessLine(".L atmFitPars.cxx+");
gROOT->ProcessLine(".L histoManager.cxx+");
gROOT->ProcessLine(".L modHistoArrayFV.cxx+");
gROOT->ProcessLine(".L mcmcApply.cxx+");
gROOT->ProcessLine(".L TH2FV.cxx+");
gROOT->ProcessLine(".L SKError.cxx+");
gROOT->ProcessLine(".L toyMC.cxx+");

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// parameters //
TString T2KMCFiles = "/nfs/data41/t2k/amissert/t2kmc/processed/feb19full_logrc/*.root";
TString MCMCFiles = "/nfs/data41/t2k/amissert/atmos/head/atmFitTools/run/results/wetrun_logrc/mcmc_run2_summary.root";
int   NMCEvents = 1e9;
//int   NMCEvents = 30000;
int   NMCMCPoints = 499;
int   index_of_pidpar = 0;
int   index_of_pi0par = 1;
int   index_of_pippar = 2;
int   index_of_rcpar  = 3;
// selection code:
//   1 -> Nu E
//   2 -> Nu Mu
//   3 -> Nu E CC1pi
int   selection_code = 3;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// run toyMC
toyMC* toy = new toyMC();
toy->indexPIDPar = index_of_pidpar;
toy->indexPi0Par = index_of_pi0par;
toy->indexPiPPar = index_of_pippar;
toy->indexRCPar  = index_of_rcpar;

// mc files
TChain* mcfiles = new TChain("h1");
mcfiles->Add(T2KMCFiles.Data());

// mcmc pars
TChain* parfiles = new TChain("MCMCpath");
parfiles->Add(MCMCFiles.Data());

toy->setChains(mcfiles,parfiles,NMCEvents);
toy->setAtmFitPars("wetrun.dat");

mcfiles->GetEntry(500); // initialize some parameters
parfiles->GetEntry(500); // initialize some parameters
//toy->makeFVUncMap(NMCMCPoints,selection_code,"./data/FVUncMapNuE1RPi.root",1);
toy->makeFVUncMap(NMCMCPoints,selection_code,"./data/FVUncMapNuE1RPiBinning2.root",2);
//toy->makeFVUncMap(NMCMCPoints,selection_code,"./data/FVUncMapNuE1RPiBinning3.root",3);
//toy->makeFVUncMap(NMCMCPoints,selection_code,"./data/FVUncMapNuEBinning2.root",2);
//toy->makeFVUncMap(NMCMCPoints,selection_code,"./data/FVUncMapNuEBinning3.root",3);

}
