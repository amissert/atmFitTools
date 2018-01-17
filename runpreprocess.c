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
 preproc->setParFileName("wetrun_newrc.dat");
 preproc->fakeShiftFlg = 0.;
 preproc->fakeNormFlg = 0.;
 preproc->runPreProcessing();

}
