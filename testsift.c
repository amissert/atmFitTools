{
 gROOT->ProcessLine(".L sift.C+");
 TChain ch("h1");
 ch.Add("/nfs/data40/t2k/amissert/skdata/atmosMC/*.reduc.09*.root");

 sift* ss = new sift(&ch);
 ss->siftIt("test"); 
}