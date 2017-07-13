{
  TString fnames = "/nfs/data41/t2k/amissert/processed/wetrun/jun27newrc/*run4_ppmc*.root";
  gROOT->ProcessLine(".L makeRCMap.C++");
  makeRCMap* makemap = new makeRCMap(fnames.Data());
  makemap->makeMap();
}

