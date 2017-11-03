{

 ///////////////////////////////////////
 // parameter file
  TString parfile = "wetrun_lowevis.dat";
  ///////////////////////////////////////

 ///////////////////////////////////////
 //load class
 gROOT->ProcessLine(".L preProcess.cxx+");

 ///////////////////////////////////////
 //setup and run preprocessing object
 preProcess* preproc = new preProcess();
// preproc->setParFileName("shimpars.dat");
// preproc->setParFileName("dryrun.dat");
 preproc->setParFileName("wetrun_lowevis.dat");
// preproc->setParFileName("allevis.dat");
// preproc->setParFileName("cosmicpars.dat");
// preproc->setRFGDir("/nfs/data41/t2k/amissert/skdata/atmospheric/sf_rfg/");
// preproc->setParFileName("t2kpars.dat");
// preproc->setParFileName("cosmicpars.dat");
// preproc->setParFileName("hpi0pars.dat");
// preproc->setParFileName("atmparsE.dat");
 preproc->fakeShiftFlg = 0.;
 preproc->fakeNormFlg = 0.;
 preproc->runPreProcessing();


}
