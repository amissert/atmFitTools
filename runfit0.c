{
  
 ///////////////////////////////////////
 // parameter file
 TString parfile = "wetrun_newrc.dat";
// TString parfile = "wetrun_lowevis.dat";
 ///////////////////////////////////////

 gROOT->ProcessLine(".L histoManager.cxx+");
 gROOT->ProcessLine(".L hSplines.cxx+");
 gROOT->ProcessLine(".L histoCompare.cxx++");
 gROOT->ProcessLine(".L atmFitPars.cxx+");
 gROOT->ProcessLine(".x ~/style.c");
 gStyle->SetOptTitle(1);
 gStyle->SetMarkerStyle(1);

  histoCompare* hc= new histoCompare(parfile.Data());
  hc->diffChainFileName = "/nfs/data41/t2k/amissert/atmos/head/atmFitTools/run/results/wetrun_newrc/atmfit*demcmc*diff*.root";
}


