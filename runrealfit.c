{
 gROOT->ProcessLine(".L histoManager.cxx+");
 gROOT->ProcessLine(".L hSplines.cxx+");
 gROOT->ProcessLine(".L histoCompare.cxx+");
 gROOT->ProcessLine(".L atmFitPars.cxx+");
 gROOT->ProcessLine(".x ~/style.c");

//  histoCompare* hc= new histoCompare("dryrun.dat");
//  histoCompare* hc= new histoCompare("wetrunrc.dat");
//  histoCompare* hc= new histoCompare("wetrun.dat");
  histoCompare* hc= new histoCompare("wetrunrc.dat");
//  hc->LnLFit();

 hc->saveFitPars("./pars/atmfit_fullpars_run3_LnLFit.root");
  



//  hc->thePars->readPars("/nfs/data41/t2k/amissert/atmos/head/atmFitTools/pars/lnlfit_priors.root");
  

//  hc->diffChainFileName = "./run/results/mcmc_diff*.root";

// hc->calcRoughParErr();
// hc->LnLFit();
}


