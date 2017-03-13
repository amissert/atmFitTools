{

gStyle->SetPalette(kLightTemperature);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// parameters //
//TString T2KMCFiles = "/Users/andy/t2k/t2kmc/processed/feb19full/*.root";
//TString T2KMCFiles = "/Users/andy/t2k/skdata/atmospheric/processed/wetrun_final/*.root";
TString T2KMCFiles = "/Users/andy/t2k/skdata/atmospheric/processed/wetrun_allevis/*.root";
TString MCMCFiles = "./results/demcmc_run2_summary.root ";
int   NMCEvents = 1e9;
//int   NMCEvents = 200000;

int   NMCMCPoints = 10;
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
//toy->modifier->flgApplyXSecPar = true;
//toy->modifier->flgApplyNormPar = false;
//toy->modifier->flgApplyFluxPar = false;
//toy->modifier->setAttFlgs(0,false);
//toy->modifier->setAttFlgs(1,false);
//toy->modifier->setAttFlgs(2,false);
//toy->modifier->setAttFlgs(3,false);
toy->fillSKErrors(NMCMCPoints,1);

}
