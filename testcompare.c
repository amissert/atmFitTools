{
 gROOT->ProcessLine(".L histoManager.C++");
 gROOT->ProcessLine(".L hSplines.C++");
 gROOT->ProcessLine(".L histoCompare.C++");
 gROOT->ProcessLine(".L atmFitPars.C++");
 gROOT->ProcessLine(".x ~/style.c");

 histoCompare* hc = new histoCompare("comptest");
 histoManager* hm = new histoManager(1000,100);
 hc->initialize(hm,hm->fitPars);
// hc->readFromFile("histos_test2",3,3,7,1);  //reads in histogram manager
// hc->readFromFile("./rootfiles/nom2_factoryOutput.root",3,3,7,1);
// hc->readFromFile(".root",3,3,7,1);

// hc->setupPars(1); //setup parameters
// hc->thePars->setSysParUnc(0,0.05);
// hc->setCompName(0,"CC1e");
// hc->setCompName(1,"CC1mu");
// hc->setCompName(2,"CCeOth");
// hc->setCompName(3,"CCmuOth");
// hc->setCompName(4,"CCOth");
// hc->setCompName(5,"NCpi0");
// hc->setCompName(6,"NCOth");
// hc->setBinName(0,"FV0");
// hc->setBinName(1,"FV1");
 //hc->setBinName(2,"FV2");
// hc->setAttName(0,"emuPID");
// hc->setAttName(0,"Other");
// hc->setRebinFactor(1);
// hc->readFitPars("./rootfiles/fitpars.root");
// hc->LnLFit();
// hc->saveFitPars("./rootfiles/fitpars_smooth.root");
// hc->addHistogram(hpidmc,0);
// hc->addHistogram(hpidmcmu,0);
// hc->addHistogram(hpiddat,1);
// hc->addHistogram(hpidmcmu,0);
// hc->addHistogram(hpiddatmu,1);

}