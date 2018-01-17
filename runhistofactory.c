{

  ///////////////////////////////////////
  // parameter file
  TString parfile = "wetrun_lowevis.dat";
  ///////////////////////////////////////
 
  //load classes
  gROOT->ProcessLine(".L histoFactory.cxx+");
  gROOT->ProcessLine(".L splineFactory.cxx+");

  // run factory
  histoFactory* hfact = new histoFactory(parfile.Data());
  hfact->runHistoFactory();

};
