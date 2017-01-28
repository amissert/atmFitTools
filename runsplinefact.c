{

  ///////////////////////////////////////
  //load classes
 // gROOT->ProcessLine(".L histoManager.cxx++");
  gROOT->ProcessLine(".L histoFactory.cxx+");
  gROOT->ProcessLine(".L splineFactory.cxx+");
  gROOT->ProcessLine(".x ~/style.c");

  //run spline factory form parameter file
  splineFactory *sfact = new splineFactory("dryrun.dat"); //< atmospheric pars
  sfact->runSplineFactory();


};
