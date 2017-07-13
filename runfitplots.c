{

 gROOT->ProcessLine(".L fitPlots.cxx++");
// gROOT->ProcessLine(".L histoManager.cxx+");
// gROOT->ProcessLine(".L hSplines.cxx+");
 gROOT->ProcessLine(".L histoCompare.cxx+");
 gROOT->ProcessLine(".L atmFitPars.cxx+");
 gROOT->ProcessLine(".L hArray.cxx+");
// gROOT->ProcessLine(".x ~/style.c");
 gStyle->SetOptTitle(1);
 gStyle->SetTitleW(0.8);
 gStyle->SetTitleH(0.06);

 histoCompare* hc= new histoCompare("wetrun_newrc.dat");

 TChain ch("MCMCpath");
// ch.Add("./run/results/wetrun/demcmc_summary_sample.root");
// ch.Add("./run/results/wetrun_logrc/mcmc*summ*.root");
 ch.Add("./run/results/wetrun_newrc/demcmc_summary_large.root");
// ch.Add("./run/results/wetrun/mcmcfit_summary_sample.root");
// ch.Add("./run/results/wetrun/demcmc_run4_summary_sample.root");

 fitPlots* fp = new fitPlots(hc,(TTree*)&ch);
 fp->nPoints = 70;
 fp->initArrays();
 fp->printFitSummaryFV("~/transfer/","demcmc_fitresult_newrc2");
// fp->fillArrays(0,0);
//fp->drawFitSummaryFV(0,0);
}
