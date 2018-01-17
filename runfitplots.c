{

 gROOT->ProcessLine(".L fitPlots.cxx++");
 gROOT->ProcessLine(".L histoCompare.cxx+");
 gROOT->ProcessLine(".L atmFitPars.cxx+");
 gROOT->ProcessLine(".L hArray.cxx+");
 gStyle->SetOptTitle(1);
 gStyle->SetTitleW(0.8);
 gStyle->SetTitleH(0.06);

 histoCompare* hc= new histoCompare("wetrun_newrc.dat");

 TChain ch("MCMCpath");
 ch.Add("./run/results/wetrun_newrc/demcmc_summary_large.root");

 fitPlots* fp = new fitPlots(hc,(TTree*)&ch);
 fp->nPoints = 70;
 fp->initArrays();
 fp->printFitSummaryFV("~/transfer/","demcmc_fitresult_newrc2");
}
