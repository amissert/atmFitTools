{

gROOT->ProcessLine(".L RCMap.C+");
TString fname = "./data/RCMap.root";
RCMap* rcmap = new RCMap(fname.Data());

}

