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

//gStyle->SetPalette(kLightTemperature);

////////////////////////////////////////////////////////////////////////////////////////////////////////
// parameters //
//TString T2KMCFiles = "/Users/andy/t2k/skdata/atmospheric/processed/wetrun_final/*.root";
TString T2KMCFiles = "/nfs/data41/t2k/amissert/processed/wetrun/allevis/*.root";
TString MCMCFiles = "./run/results/wetrun_logrc/demcmc_run2_summary.root";
int   NMCEvents = 1e9;
//int   NMCEvents = 200000;
int   NMCMCPoints = 5000;
int   NMarginalPoints = 10;
int   index_of_pidpar = 0;
int   index_of_pi0par = 1;
int   index_of_pippar = 2;
int   index_of_pmom   = -1;
int   index_of_rcpar  = 3;
////////////////////////////////////////////////////////////////////////////////////////////////////////


// run toyMC
toyMC* toy = new toyMC();
toy->indexPIDPar = index_of_pidpar;
toy->indexPi0Par = index_of_pi0par;
toy->indexPiPPar = index_of_pippar;
toy->indexRCPar  = index_of_rcpar;
//toy->indexMom    = index_of_pmom;

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

//toy->fillMarginalizedSKErr(NMCMCPoints,NMarginalPoints);

toy->modifier->flgGlobalUseBestPars = true;
toy->modifier->setUseBestFitSystPars(true);
toy->modifier->setUseBestFitNormPars(true);
//toy->modifier->setAllAlphaFlags(false);
//toy->skErr->effDefinition = 1;
//toy->fillSKErrors(NMCMCPoints,1,0,0);
toy->fillSKErrors(NMCMCPoints,1,1,0);

//toy->skErr->drawVariable(toy->fastevents,"fqpidpar",1,1)
}
