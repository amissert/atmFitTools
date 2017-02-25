{
gROOT->ProcessLine(".L histoCompare.cxx+");
gROOT->ProcessLine(".L atmFitPars.cxx+");
gROOT->ProcessLine(".L histoManager.cxx+");
gROOT->ProcessLine(".L modHistoArray.cxx+");
gROOT->ProcessLine(".L modHistoArrayFV.cxx+");
gROOT->ProcessLine(".L mcmcApply.cxx+");
gROOT->ProcessLine(".L TH2FV.cxx+");
gROOT->ProcessLine(".L SKError.cxx+");
gROOT->ProcessLine(".L toyMC.cxx+");

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// parameters //
//TString T2KMCFiles = "/nfs/data41/t2k/amissert/t2kmc/processed/jan30full//*.root";
//TString T2KMCFiles = "/nfs/data41/t2k/amissert/t2kmc/processed/feb6full/*.root";
TString T2KMCFiles = "/nfs/data41/t2k/amissert/t2kmc/processed/feb19full_logrc/*.root";
//TString MCMCFiles = "/nfs/data41/t2k/amissert/atmos/head/atmFitTools/run/results/wetrun/demcmc_summary_sample.root";
//TString MCMCFiles = "/nfs/data41/t2k/amissert/atmos/head/atmFitTools/run/results/wetrun/demcmc_run5_summary_sample.root";
TString MCMCFiles = "/nfs/data41/t2k/amissert/atmos/head/atmFitTools/run/results/wetrun_logrc/mcmc_run2_summary.root";
//TString MCMCFiles = "/nfs/data41/t2k/amissert/atmos/head/atmFitTools/run/results/wetrun/mcmcfit_run1_summary_sample.root";
//int   NMCEvents = 1e9;
int   NMCEvents = 200000;
int   NMCMCPoints = 100;
int   index_of_pidpar = 0;
int   index_of_pi0par = 1;
int   index_of_pippar = 2;
int   index_of_pmom   = -1;
int   index_of_rcpar  = 3;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// run toyMC
toyMC* toy = new toyMC();
toy->indexPIDPar = index_of_pidpar;
toy->indexPi0Par = index_of_pi0par;
toy->indexPiPPar = index_of_pippar;
toy->indexRCPar  = index_of_rcpar;
toy->indexMom    = index_of_pmom;

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
toy->makeFVMapNuMu(NMCMCPoints,"./data/FVUncMapNuMu.root");
//toy->makeFVMapNuE(NMCMCPoints,"./data/FVUncMapNuE.root");
//toy->makeFVMapNuE(NMCMCPoints);
//toy->hArrFV->saveClose();

toy->fillSKErrors(100);

}
