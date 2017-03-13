{
  TH1D* hh = new TH1D("hh","hh",10,0,10);
  hh->Fill(7);
  hh->Fill(7);
  hArray* ha = new hArray("tst",hh,10);
}
