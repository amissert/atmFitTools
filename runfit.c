{
 gROOT->ProcessLine(".L histoManager.cxx+");
 gROOT->ProcessLine(".L hSplines.cxx+");
 gROOT->ProcessLine(".L histoCompare.cxx++");
 gROOT->ProcessLine(".L atmFitPars.cxx+");
 gROOT->ProcessLine(".x ~/style.c");
 gStyle->SetOptTitle(1);
 gStyle->SetMarkerStyle(1);

// Initialize with parameter file name
histoCompare* hc= new histoCompare("<parfile>.dat");


// Can optionally fix some parameters in the fit by changing the fit flags...
// hc->thePars->fixAllAttPars(0);
// hc->thePars->fixAllAttPars(1);
// hc->thePars->fixAllAttPars(2);

// Can optionally set physical boundaries for some parameters (not needed with new RC definition)
//  hc->hManager->setLoBound(3,0);

// You can also load previously saved parameter  
//    hc->thePars->readPars("./run/results/demcmc_testmcmcfit_41.root");

// Option to change the spline flag (now in parameters file).
// hc->hManager->useSplineFlg = 0;

// This runs a fit using minuit
// hc->LnLFit();

// Saves pars so they can be loaded later
// hc->saveFitPars("./pars/fullpars_wetrun_LnLfit_logrcpar.root");

// Set the tuning parameter that determines mcmc step size
// hc->tunePar = 0.033;

// run Markov chain
// hc->runMCMC(2000000);

// For DE-MCMC, specify file with differential steps
// hc->diffChainFileName = "mcmcdiff.root";
// hc->runDEMCMC(2000000);

}


