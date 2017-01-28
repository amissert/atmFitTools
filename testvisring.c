{
 gROOT->ProcessLine(".L visRing.cxx+");
 gROOT->ProcessLine(".L fqEvent.cxx+");


 TChain ch("h1");
 ch.Add("/nfs/data41/t2k/amissert/skdata/atmospheric/shimp_mc/mar16sk4.reduc.0488_fQv5r0.root");
 
 fqEvent* evt = new fqEvent(&ch);

 visRing* vis = new visRing(evt);
 evt->GetEntry(158);

}
