{
 gROOT->ProcessLine(".L histoManager.cxx+");
 gROOT->ProcessLine(".L hSplines.cxx+");
 gROOT->ProcessLine(".L histoCompare.cxx++");
 gROOT->ProcessLine(".L atmFitPars.cxx+");
 gROOT->ProcessLine(".x ~/style.c");
 gStyle->SetOptTitle(1);
 gStyle->SetMarkerStyle(1);

//  histoCompare* hc= new histoCompare("dryrun.dat");
//  histoCompare* hc= new histoCompare("wetrun.dat");
  histoCompare* hc= new histoCompare("wetrun_2bin.dat");

// hc->thePars->fixAllAttPars(0);
// hc->thePars->fixAllAttPars(1);
// hc->thePars->fixAllAttPars(2);

//  hc->hManager->setLoBound(3,0);


//  histoCompare* hc= new histoCompare("shimpars.dat");

//  hc->MCMCOutputFile = outfilename.Data();
//  hc->hManager->setLoBound(3,0);
//  hc->hManager->setLoBound(2,0);
//  hc->thePars->readPars("/nfs/data41/t2k/amissert/atmos/head/atmFitTools/pars/lnlfit_priors.root");
  
//    hc->thePars->readPars("./run/results/demcmc_testmcmcfit_41.root");

//  hc->thePars->readPars("/nfs/data41/t2k/amissert/atmos/head/atmFitTools/pars/initial_lnlfit_pars.root");
//  hc->setPar(26,1.0);
//  hc->setPar(27,0.0); 
//  hc->setPar(154,1.0);
//  hc->setPar(155,0.0);
//  hc->setPar(106,1.0);
//  hc->setPar(107,0.0);
//  hc->thePars->fixPar[26]=1;
//  hc->thePars->fixPar[27]=1;
//  hc->thePars->fixPar[154]=1;
//  hc->thePars->fixPar[155]=1;
//  hc->thePars->fixPar[106]=1;
//  hc->thePars->fixPar[107]=1;
//  hc->tunePar = 0.02;

//  hc->diffChainFileName = "./run/results/mcmc_diff*.root";

// histoCompare* hc= new histoCompare("atmparsE.dat");
// histoCompare* hc= new histoCompare("fakepars2.dat");
// histoCompare* hc= new histoCompare("fakepars1.dat");
// histoCompare* hc= new histoCompare("shimpars.dat");
// hc->hManager->setLoBound(2,0);
// hc->hManager->setLoBound(3,0);
// hc->tunePar = 0.01;
// hc->thePars->readPars("./pars/inital_lnlfit_pars.root");
// hc->LnLFit();
// hc->thePars->readPars("/nfs/data41/t2k/amissert/atmos/head/atmFitTools/initial_lnlfit.root");
// hc->thePars->readPars("/nfs/data41/t2k/amissert/atmos/head/atmFitTools/run/results/demcmc_testmcmcfit_2.root");

// hc->thePars->fixAllSystPars(1);

// hc->hManager->useSplineFlg = 0;
// hc->LnLFit();
// hc->saveFitPars("./pars/fullpars_wetrun_LnLfit_logrcpar.root");
// hc->tunePar = 0.033;
// hc->runMCMC(2000000);

// hc->thePars->readPars("./pars/fullfit_nosyst.root");


// hc->diffChainFileName = "mcmcdiff.root";

 //hc->thePars->fixAllAttPars(2);
 //hc->thePars->fixAllAttPars(3);
 //hc->thePars->fixAllAttPars(4);
// hc->calcRoughParErr();
// hc->LnLFit();
}


