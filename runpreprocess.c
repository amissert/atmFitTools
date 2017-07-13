{

 ///////////////////////////////////////
 // parameter file
 TString parfile = "wetrun_lowevis.dat";
 ///////////////////////////////////////

 //load class
 gROOT->ProcessLine(".L preProcess.cxx+");

 //setup and run preprocessing
 preprocsetParFileName(parfile.Data());
 preproc->fakeShiftFlg = 0.;
 preproc->fakeNormFlg = 0.;
 preproc->runPreProcessing();


}
