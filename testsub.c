{


//gStyle->SetPalette(kTemperatureMap);
gROOT->ProcessLine(".x ~/corr.c");
 gROOT->ProcessLine(".L makeCov.cxx+");
 gStyle->SetOptStat(0);
 gStyle->SetOptTitle(0);
// TFile f("./mcmc/closure_allpar_mcmcfit.root");
// TFile f("./run/results/wetrun/demcmc_.root");
// TFile f("./run/results/wetrun/demcmc_summary_sample.root");
// TFile f("./run/results/wetrun/demcmc_run5_summary_sample.root");
 TFile f("/nfs/data41/t2k/amissert/atmos/head/atmFitTools/run/results/wetrun_logrc/demcmc_large_summary.root");
// TFile f("./run/results/wetrun/mcmcfit_run5_summary_sample.root");
 TTree* tr = (TTree*)f.Get("MCMCpath");
 makeCov *maker = new makeCov("wetrun.dat");
// makeCov *maker = new makeCov("fakepars1.dat");
 maker->setParTree(tr);
 maker->nburn = 1000;
 maker->nskip = 50;
 //gStyle->SetPalette(kBlackBody);
 maker->buildMatrix();
 maker->drawLabeldCor();
 maker->drawBinSubMatrix(2);
 maker->drawSysValues();
 maker->drawBinSubMatrix(5,1);
}
