     
{
  // corr pallette
  TCanvas *c2  = new TCanvas("c2","c2",0,0,600,400);
  TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",1,3,1,3);
  const Int_t Number = 5;
  Double_t Red[Number]    = {0.00, 0.00, 1.00, 1.00, 1.00};
  Double_t Green[Number]  = {0.00, 1.00, 1.00, 1.00, 0.00};
  Double_t Blue[Number]   = {1.00, 1.00, 1.00, 0.00, 0.00};
  Double_t Length[Number] = {0.00, 0.40, 0.50, 0.60, 1.00};
  Int_t nb=50;
  TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
  f2->SetContour(nb);
  f2->SetLineWidth(1);
  f2->SetLineColor(kBlack);
  f2->Draw("surf1z");
  return c2;
}
