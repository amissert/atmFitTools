{

  ///////////////////////////////////////
  // parameter file
  TString parfile = "wetrun_lowevis.dat";
  ///////////////////////////////////////

  //load classes
  gROOT->ProcessLine(".L histoFactory.cxx+");
  gROOT->ProcessLine(".L splineFactory.cxx+");
  gROOT->ProcessLine(".x ~/style.c");

  //run spline factory form parameter file
  splineFactory *sfact = new splineFactory(parfile.Data()); //< atmospheric pars
  sfact->runSplineFactory();

};
