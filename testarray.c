{
gROOT->ProcessLine(".L modHistoArrayFV.cxx+");
gROOT->ProcessLine(".L TH2FV.cxx+");
modHistoArrayFV* mm = new modHistoArrayFV("NuMu_histogram_array.root",6,100);
//modHistoArrayFV* mm = new modHistoArrayFV("hE_nuE_v2_histogram_array.root",6,100);
mm->calcSummary();
//mm->printUncMap("~/transfer/");
}
