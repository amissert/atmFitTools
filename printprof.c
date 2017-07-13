{
 gROOT->ProcessLine(".L histoManager.cxx+");
 gROOT->ProcessLine(".L hSplines.cxx+");
 gROOT->ProcessLine(".L histoCompare.cxx++");
 gROOT->ProcessLine(".L atmFitPars.cxx+");
 gROOT->ProcessLine(".x ~/style.c");
 gStyle->SetOptTitle(1);
 gStyle->SetMarkerStyle(1);

 /////////////////////////////////////////////////
 double xrange = 100;
// TString parfile = "wetrun.dat";
// TString outfile = "~/prof_many_mu_mupipid_dr5.root";
 TString xname = "#mu/#pi bias for single e";
 TString parfile = "wetrun_2bin.dat";
 TString outfile = "~/prof_4bin_mu_mupipid_dr5.root";
 int icomp=1;
 int ibin=5;
 int iatt=2;
 int itype=1;
 /////////////////////////////////////////////////
 histoCompare* hc= new histoCompare(parfile.Data());
 histoCompare* hc2= new histoCompare(parfile2.Data());
 int index = hc->thePars->getParIndex(ibin,icomp,iatt,itype);
 hc->thePars->fixAllPars(1);
// hc->thePars->fixPars[index]=0;
 hc->thePars->fixPar[index] = 0;
// hc->LnLFit();
 hc->profileL(index, xrange);
 hc->hProf[0]->GetXaxis()->SetTitle(xname.Data());
 hc->hProf[0]->SaveAs(outfile.Data());
// hc2->thePars->fixAllPars(1);
// hc2->thePars->fixPars[index]=0;
// hc2->LnLFit();
// hc2->profileL(index, xrange);
// hc2->hProf[0]->SaveAs(outfile.Data());
// hc->hProf[0]->Draw("");
// hc2->hProf[0]->Draw("same");

}


