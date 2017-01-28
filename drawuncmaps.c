{

TString pltname = plotdir.Data();
pltname.Append(hccqe->GetName());
pltname.Append(".png");
hccqe->Draw("colz");
c1->Print(pltname.Data());
//
TString pltname = plotdir.Data();
pltname.Append(hccnqe->GetName());
pltname.Append(".png");
hccnqe->Draw("colz");
c1->Print(pltname.Data());
//
TString pltname = plotdir.Data();
pltname.Append(hccwrong->GetName());
pltname.Append(".png");
hccwrong->Draw("colz");
c1->Print(pltname.Data());
//
TString pltname = plotdir.Data();
pltname.Append(hnc->GetName());
pltname.Append(".png");
hnc->Draw("colz");
c1->Print(pltname.Data());
//




}
