{

//
TFile f("./data/FVUncMapNuE.root");
//TFile f("./data/FVUncMapNuMu.root");

// read in histograms
TH2Poly* hccqe =  (TH2Poly*)f.Get("FVUncMapCCQE");
TH2Poly* hccnqe = (TH2Poly*)f.Get("FVUncMapCCnQE");
TH2Poly* hccwrong = (TH2Poly*)f.Get("FVUncMapCCWrong");
TH2Poly* hnc = (TH2Poly*)f.Get("FVUncMapNC");
TString plotdir = "~/transfer/";

hccqe->GetXaxis()->SetRangeUser(0,1000);
hccnqe->GetXaxis()->SetRangeUser(0,1000);
hccwrong->GetXaxis()->SetRangeUser(0,1000);
hnc->GetXaxis()->SetRangeUser(0,1000);

hccqe->GetYaxis()->SetRangeUser(0,1200);
hccnqe->GetYaxis()->SetRangeUser(0,1200);
hccwrong->GetYaxis()->SetRangeUser(0,1200);
hnc->GetYaxis()->SetRangeUser(0,1200);

hccqe->SetStats(0);
hccnqe->SetStats(0);
hccwrong->SetStats(0);
hnc->SetStats(0);
hccqe->SetTitle(0);
hccnqe->SetTitle(0);
hccwrong->SetTitle(0);
hnc->SetTitle(0);

gROOT->ProcessLine(".x drawuncmaps.c");

}
