{
gStyle->SetPalette(kTemperatureMap);
// gROOT->ProcessLine(".L makeCov.cxx+");
// TFile f("./mcmc/closure_allpar_mcmcfit.root");
// TFile f("./run/results/wetrun/demcmc_.root");
// TFile f("./run/results/wetrun/demcmc_summary_sample.root");
// TFile f("./run/results/wetrun/demcmc_run5_summary_sample.root");
 TFile f("./results/demcmc_run2_summary.root");
// TFile f("./run/results/wetrun/mcmcfit_run5_summary_sample.root");
 TTree* tr = (TTree*)f.Get("MCMCpath");
 makeCov *maker = new makeCov("wetrun.dat");
// makeCov *maker = new makeCov("fakepars1.dat");
 maker->setParTree(tr);
 maker->nburn = 400;
 //gStyle->SetPalette(kBlackBody);
 maker->buildMatrix();
}
