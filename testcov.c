{

 ///////////////////////////////////////
 // parameter file
 TString parfile = "wetrun_newrc.dat";
 TString mcmcparfile = "/nfs/data41/t2k/amissert/atmos/head/atmFitTools/run/results/wetrun_newrc/demcmc_summary_run9.root";
// TString mcmcparfile = "/nfs/data41/t2k/amissert/atmos/head/atmFitTools/run/results/wetrun_lowevis/mcmc_summary.root";
// TString mcmcparfile = "./run/results/wetrun_logrc/demcmc_large_summary.root";
 ///////////////////////////////////////


 // make covariance matricies
 gROOT->ProcessLine(".x ~/corr.c");
 gROOT->ProcessLine(".L makeCov.cxx++");
 TFile f(mcmcparfile.Data());
 TTree* tr = (TTree*)f.Get("MCMCpath");
 makeCov *maker = new makeCov(parfile.Data());
 maker->setParTree(tr);
 maker->nburn = 50;
 maker->nskip = 1;
 maker->buildMatrix();
 maker->drawLabeldCor();

}
