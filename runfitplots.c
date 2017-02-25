{

 gROOT->ProcessLine(".L fitPlots.cxx++");
// gROOT->ProcessLine(".L histoManager.cxx+");
// gROOT->ProcessLine(".L hSplines.cxx+");
 gROOT->ProcessLine(".L histoCompare.cxx++");
 gROOT->ProcessLine(".L atmFitPars.cxx++");
 gROOT->ProcessLine(".x ~/style.c");

 histoCompare* hc= new histoCompare("wetrun.dat");

 TChain ch("MCMCpath");
// ch.Add("./run/results/wetrun/demcmc_summary_sample.root");
 ch.Add("./run/results/wetrun_logrc/mcmc*summ*.root");
// ch.Add("./run/results/wetrun/mcmcfit_summary_sample.root");
// ch.Add("./run/results/wetrun/demcmc_run4_summary_sample.root");

 fitPlots* fp = new fitPlots(hc,(TTree*)&ch);
 fp->nPoints = 100;
 fp->initArrays();
 fp->printFitSummary("~/transfer/","wetfit");
// fp->fillArrays(0,0);

}
