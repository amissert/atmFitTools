{
gROOT->ProcessLine(".L moreUncertainties.cxx+");
gROOT->ProcessLine(".L fqProcessedEvent.cxx+");
gROOT->ProcessLine(".L TH2FV.cxx+");


 moreUncertainties* more = new moreUncertainties("/nfs/data41/t2k/amissert/atmos/head/atmFitTools/data/");

 TChain ch("h1");
 ch.Add("/nfs/data41/t2k/amissert/t2kmc/processed/nov29/*8.root");
 fqProcessedEvent* fq = new fqProcessedEvent(&ch);

 more->setEventPointer(fq);

// more->fillFVHisto(&ch);
 
}
