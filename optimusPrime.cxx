#ifndef PRIME_CXX
#define PRIME_CXX


#include "optimusPrime.h"



////////////////////////////////////////////////////////////////////
void optimusPrime::drawFOM(int nselection, int noscpar){

  TString sampname[] = {"#nu_{e} CCQE","#nu_{#mu} CCQE", "#nu_{e} CC1R#pi"};
  TString osparname[] = {"#theta_{23}","#Delta M^{2}","#theta_{13}","#delta_{CP}"};
 // 0 -> figure of merit
 // 1 -> dN/dtheta
 // 2 -> (dN/dtheta) / N
 // 3 -> N
 // 4 -> sys / N
  TString ztitle[] = {"Figure of Merit","|dN/d#theta|","Avg. |dN/d#theta|",
                      "N","Avg. #sigma_{sys}"};
  TCanvas* cc = new TCanvas("cc","cc",800,700);
  cc->GetPad(0)->SetLeftMargin(0.10);
  cc->GetPad(0)->SetRightMargin(0.18);
  cc->GetPad(0)->SetBottomMargin(0.10);
  hFV->GetXaxis()->SetTitle("Towall Cut [cm]");
  hFV->GetYaxis()->SetTitle("Wall Cut [cm]");
  hFV->GetZaxis()->SetTitle(ztitle[FOMType]);
  hFV->GetZaxis()->SetTitleOffset(1.5);
  hFV->GetZaxis()->SetNdivisions(5);
  hFV->SetTitle(sampname[nselection-1].Data());
  smartAxis(hFV);
  hFV->Draw("colz");
  hFV->Draw("colz");

  return;
}



////////////////////////////////////////////////////////////////////
void optimusPrime::setSeedHisogram(int nselection, TH1D* histo, int oscpar){

  if (nselection==2){
    if (oscpar==0){
      hErecSeedNuMuTh23 = histo;
    }
    if (oscpar==1){
      hErecSeedNuMudM2 = histo;
    }
  }
  if (nselection==1){
    hErecSeedNuE = histo;
  }
  if (nselection==3){
    hErecSeedNuE = histo;
  }

  return;
}



////////////////////////////////////////////////////////////////////
void optimusPrime::printCompare(const char* dir,float tw1, float w1, float tw2, float w2, int oscpar, int nselection){

  // plot name setup
  TString plotdir = dir;
  TString name = "";

  // compare plots
  compareCuts(tw1, w1, tw2, w2, oscpar, nselection);

  // print it
  name = plotdir.Data();
  if (nselection==2) name.Append("numu_");
  else if (nselection==1){
    name.Append("nue_");
  }
  else if (nselection==3){
    name.Append("nuecc1pi_");
  }

  name.Append(Form("compare_%dw%d_vs_%dw%d.png",(int)tw1,(int)w1,(int)tw2,(int)w2));
  multiPad->Print(name.Data());

  // print it in words
  canPad = new TCanvas("can","can",800,700);
  printCutDiff(nselection);
  name = plotdir.Data();
  if (nselection==2) name.Append("numu_");
  else if (nselection==1){
    name.Append("nue_");
  }
  else if (nselection==3){
    name.Append("nuecc1pi_");
  }
  name.Append(Form("text_%dw%d_vs_%dw%d.png",(int)tw1,(int)w1,(int)tw2,(int)w2));
  canPad->Print(name.Data());

  return;
}



////////////////////////////////////////////////////////////////////
void optimusPrime::printCutDiff(int nselection){

  txtSummary = new TPaveText(0.1,0.1,0.9,0.9,"NB");
  txtSummary->SetFillColor(kWhite);

  const int nlines = 6;
  TString txtline[nlines];

  int iline =0;
  if (nselection==2){
    float n1;
    float n2;
    float diff;
    float pdiff;
    // get numbers
    n1 = plots1->pltEnuMu->Integral(1, plots1->pltEnuMu->GetNbinsX()+1);
    n2 = plots2->pltEnuMu->Integral(1,plots2->pltEnuMu->GetNbinsX()+1);
    diff = n2 - n1;
    if (n1>0.){
      pdiff = 100.*(n2 - n1)/n1;
    }
    else{
      pdiff = 0.;
    }
    txtline[0] = Form("Total: %5.2f #rightarrow %5.2f #Delta = %5.2f (%5.2f \%)",n1,n2,diff,pdiff);
    txtSummary->AddText(txtline[0].Data());
    iline++;
    // 
    n1 = plots1->pltEnuMuCat[0]->Integral(1,plots1->pltEnuMuCat[0]->GetNbinsX()+1);
    n2 = plots2->pltEnuMuCat[0]->Integral(1,plots2->pltEnuMuCat[0]->GetNbinsX()+1);
    diff = n2 - n1;
    if (n1>0.){
      pdiff = 100.*(n2 - n1)/n1;
    }
    else{
      pdiff = 0.;
    }
    txtline[iline] = Form("CCQE: %5.2f #rightarrow %5.2f #Delta = %5.2f (%5.2f \%)",n1,n2,diff,pdiff);
    txtSummary->AddText(txtline[iline].Data());
    iline++;
    // 
    n1 = plots1->pltEnuMuCat[1]->Integral(1,plots1->pltEnuMuCat[1]->GetNbinsX()+1);
    n2 = plots2->pltEnuMuCat[1]->Integral(1,plots2->pltEnuMuCat[1]->GetNbinsX()+1);
    diff = n2 - n1;
    if (n1>0.){
      pdiff = 100.*(n2 - n1)/n1;
    }
    else{
      pdiff = 0.;
    }
    txtline[iline] = Form("CCnQE: %5.2f #rightarrow %5.2f #Delta = %5.2f (%5.2f \%)",n1,n2,diff,pdiff);
    txtSummary->AddText(txtline[iline].Data());
    iline++;
    // 
    n1 = plots1->pltEnuMuCat[2]->Integral(1, plots1->pltEnuMuCat[2]->GetNbinsX()+1);
    n2 = plots2->pltEnuMuCat[2]->Integral(1, plots1->pltEnuMuCat[2]->GetNbinsX()+1);
    diff = n2 - n1;
    if (n1>0.){
      pdiff = 100.*(n2 - n1)/n1;
    }
    else{
      pdiff = 0.;
    }
    txtline[iline] = Form("CCMisID: %5.2f #rightarrow %5.2f #Delta = %5.2f (%5.2f \%)",n1,n2,diff,pdiff);
    txtSummary->AddText(txtline[iline].Data());
    iline++;
    // 
    n1 = plots1->pltEnuMuCat[3]->Integral(1,plots1->pltEnuMuCat[3]->GetNbinsX()+1);
    n2 = plots2->pltEnuMuCat[3]->Integral(1,plots2->pltEnuMuCat[3]->GetNbinsX()+1);
    diff = n2 - n1;
    if (n1>0.){
      pdiff = 100.*(n2 - n1)/n1;
    }
    else{
      pdiff = 0.;
    }
    txtline[iline] = Form("NC: %5.2f #rightarrow %5.2f #Delta = %5.2f (%5.2f \%)",n1,n2,diff,pdiff);
    txtSummary->AddText(txtline[iline].Data());
    iline++;
    // 
    n1 = plots1->pltEnuMuCat[4]->Integral(1,plots2->pltEnuMuCat[4]->GetNbinsX()+1);
    n2 = plots2->pltEnuMuCat[4]->Integral(1,plots2->pltEnuMuCat[4]->GetNbinsX()+1);
    diff = n2 - n1;
    if (n1>0.){
      pdiff = 100.*(n2 - n1)/n1;
    }
    else{
      pdiff = 0.;
    }
    txtline[iline] = Form("Entering: %5.2f #rightarrow %5.2f #Delta = %5.2f (%5.2f \%)",n1,n2,diff,pdiff);
    txtSummary->AddText(txtline[iline].Data());
    iline++;
  }

 else if (nselection==1) {
    float n1;
    float n2;
    float diff;
    float pdiff;
    // get numbers
    n1 = plots1->pltEnuE->Integral();
    n2 = plots2->pltEnuE->Integral();
    diff = n2 - n1;
    if (n1>0.){
      pdiff = 100.*(n2 - n1)/n1;
    }
    else{
      pdiff = 0.;
    }
    txtline[0] = Form("Total: %5.2f #rightarrow %5.2f #Delta = %5.2f (%5.2f \%)",n1,n2,diff,pdiff);
    txtSummary->AddText(txtline[0].Data());
    iline++;
    // 
    n1 = plots1->pltEnuECat[0]->Integral();
    n2 = plots2->pltEnuECat[0]->Integral();
    diff = n2 - n1;
    if (n1>0.){
      pdiff = 100.*(n2 - n1)/n1;
    }
    else{
      pdiff = 0.;
    }
    txtline[iline] = Form("CCQE: %5.2f #rightarrow %5.2f #Delta = %5.2f (%5.2f \%)",n1,n2,diff,pdiff);
    txtSummary->AddText(txtline[iline].Data());
    iline++;
    // 
    n1 = plots1->pltEnuECat[1]->Integral();
    n2 = plots2->pltEnuECat[1]->Integral();
    diff = n2 - n1;
    if (n1>0.){
      pdiff = 100.*(n2 - n1)/n1;
    }
    else{
      pdiff = 0.;
    }
    txtline[iline] = Form("CCnQE: %5.2f #rightarrow %5.2f #Delta = %5.2f (%5.2f \%)",n1,n2,diff,pdiff);
    txtSummary->AddText(txtline[iline].Data());
    iline++;
    // 
    n1 = plots1->pltEnuECat[2]->Integral();
    n2 = plots2->pltEnuECat[2]->Integral();
    diff = n2 - n1;
    if (n1>0.){
      pdiff = 100.*(n2 - n1)/n1;
    }
    else{
      pdiff = 0.;
    }
    txtline[iline] = Form("CCMisID: %5.2f #rightarrow %5.2f #Delta = %5.2f (%5.2f \%)",n1,n2,diff,pdiff);
    txtSummary->AddText(txtline[iline].Data());
    iline++;
    // 
    n1 = plots1->pltEnuECat[3]->Integral();
    n2 = plots2->pltEnuECat[3]->Integral();
    diff = n2 - n1;
    if (n1>0.){
      pdiff = 100.*(n2 - n1)/n1;
    }
    else{
      pdiff = 0.;
    }
    txtline[iline] = Form("NC: %5.2f #rightarrow %5.2f #Delta = %5.2f (%5.2f \%)",n1,n2,diff,pdiff);
    txtSummary->AddText(txtline[iline].Data());
    iline++;
    // 
    n1 = plots1->pltEnuECat[4]->Integral();
    n2 = plots2->pltEnuECat[4]->Integral();
    diff = n2 - n1;
    if (n1>0.){
      pdiff = 100.*(n2 - n1)/n1;
    }
    else{
      pdiff = 0.;
    }
    txtline[iline] = Form("Entering: %5.2f #rightarrow %5.2f #Delta = %5.2f (%5.2f \%)",n1,n2,diff,pdiff);
    txtSummary->AddText(txtline[iline].Data());
    iline++;
  }

 else if (nselection==3) {
    float n1;
    float n2;
    float diff;
    float pdiff;
    // get numbers
    n1 = plots1->pltEnuE1RPi->Integral();
    n2 = plots2->pltEnuE1RPi->Integral();
    diff = n2 - n1;
    if (n1>0.){
      pdiff = 100.*(n2 - n1)/n1;
    }
    else{
      pdiff = 0.;
    }
    txtline[0] = Form("Total: %5.2f #rightarrow %5.2f #Delta = %5.2f (%5.2f \%)",n1,n2,diff,pdiff);
    txtSummary->AddText(txtline[0].Data());
    iline++;
    // 
    n1 = plots1->pltEnuE1RPiCat[0]->Integral();
    n2 = plots2->pltEnuE1RPiCat[0]->Integral();
    diff = n2 - n1;
    if (n1>0.){
      pdiff = 100.*(n2 - n1)/n1;
    }
    else{
      pdiff = 0.;
    }
    txtline[iline] = Form("CCQE: %5.2f #rightarrow %5.2f #Delta = %5.2f (%5.2f \%)",n1,n2,diff,pdiff);
    txtSummary->AddText(txtline[iline].Data());
    iline++;
    // 
    n1 = plots1->pltEnuE1RPiCat[1]->Integral();
    n2 = plots2->pltEnuE1RPiCat[1]->Integral();
    diff = n2 - n1;
    if (n1>0.){
      pdiff = 100.*(n2 - n1)/n1;
    }
    else{
      pdiff = 0.;
    }
    txtline[iline] = Form("CCnQE: %5.2f #rightarrow %5.2f #Delta = %5.2f (%5.2f \%)",n1,n2,diff,pdiff);
    txtSummary->AddText(txtline[iline].Data());
    iline++;
    // 
    n1 = plots1->pltEnuE1RPiCat[2]->Integral();
    n2 = plots2->pltEnuE1RPiCat[2]->Integral();
    diff = n2 - n1;
    if (n1>0.){
      pdiff = 100.*(n2 - n1)/n1;
    }
    else{
      pdiff = 0.;
    }
    txtline[iline] = Form("CCMisID: %5.2f #rightarrow %5.2f #Delta = %5.2f (%5.2f \%)",n1,n2,diff,pdiff);
    txtSummary->AddText(txtline[iline].Data());
    iline++;
    // 
    n1 = plots1->pltEnuE1RPiCat[3]->Integral();
    n2 = plots2->pltEnuE1RPiCat[3]->Integral();
    diff = n2 - n1;
    if (n1>0.){
      pdiff = 100.*(n2 - n1)/n1;
    }
    else{
      pdiff = 0.;
    }
    txtline[iline] = Form("NC: %5.2f #rightarrow %5.2f #Delta = %5.2f (%5.2f \%)",n1,n2,diff,pdiff);
    txtSummary->AddText(txtline[iline].Data());
    iline++;
    // 
    n1 = plots1->pltEnuE1RPiCat[4]->Integral();
    n2 = plots2->pltEnuE1RPiCat[4]->Integral();
    diff = n2 - n1;
    if (n1>0.){
      pdiff = 100.*(n2 - n1)/n1;
    }
    else{
      pdiff = 0.;
    }
    txtline[iline] = Form("Entering: %5.2f #rightarrow %5.2f #Delta = %5.2f (%5.2f \%)",n1,n2,diff,pdiff);
    txtSummary->AddText(txtline[iline].Data());
    iline++;
  }
  txtSummary->Draw();
}




////////////////////////////////////////////////////////////////////////
// make FV map
// make Delta plots
// apply Delta plots as mask
// save all histograms
void optimusPrime::makeAllPlots(float twmax, float wmax, int oscpar, int npts, int nselection){


  // name setup///////////////////////////////////
  if (!outDir.CompareTo("")){
    outDir = "~/transfer/";
  }
  pltTag = "plot_";
  if (nselection==1){
    pltTag.Append("nuE");
  }
  else if (nselection==2){
    pltTag.Append("nuMu");
  }
  else if (nselection==3){
    pltTag.Append("nuE1Rpi");
  }


  TString prefix = outDir.Data();
  prefix.Append(pltTag.Data());
  TString pltname;
  TString pltfile;

  // canvas setup//////////////////////////////////
  TCanvas* cc = new TCanvas("cc","cc",800,700);

  // get map of figure of merit and calc best bin//
  cout<<"Calculating FOM Map: "<<endl;
  calcFOMMap(twmax,wmax,oscpar,npts,nselection);

  // print figure of merit and save ///////////////
  hFV->SetStats(0);
  hFV->SetTitle(0);
  hFV->Draw("colz");
  pltname = prefix.Data();
  pltname.Append("_FOM.png");
  // plot names
  cc->Print(pltname.Data());
  pltfile = prefix.Data();
  pltfile.Append("_FOM.root");
  hFV->SaveAs(pltfile.Data());

  // get 5% variation in systematic estimat
  cout<<"Defining small variation: "<<endl;
  double besttowall = hFV->GetBinCenterX(bestFOMbin);
  double bestwall = hFV->GetBinCenterY(bestFOMbin);
  flgPrintSummary = 1;
  calcFOMBinned(nselection,besttowall,bestwall,oscpar,1);
  cout<<"Small Varation: "<<smallVariation<<endl;
  flgPrintSummary = 0;

  // setup a mask for small variation
  hMask = (TH2FV*)hFV->Clone("hmask");
  maskThresh = bestFOMvalue - smallVariation;
  hMask->SetMinimum(maskThresh);
  hMask->SetMaximum(bestFOMvalue);
  flgUseMask = 1.;

  // now print out the mask:
  for (int i=0; i<hMask->GetNumberOfBins(); i++){
    double binc = hMask->GetBinContent(i);
    if (binc<maskThresh) hMask->SetBinContent(i,0.);
    else{
      hMask->SetBinContent(i,1.);
    }
  }
  maskThresh = 0.8;
  hMask->SetMinimum(0.9);
  hMask->SetMaximum(1.2);
  hMask->SetStats(0);
  hMask->SetTitle(0);
  hMask->Draw("col");
  pltname = prefix.Data();
  pltname.Append("_Preferred_Region.png");
  cc->Print(pltname.Data());
  pltfile = prefix.Data();
  pltfile.Append("_Preferred_Region.root");  
  hMask->SaveAs(pltfile.Data());
 
  // turn mask back off
  flgUseMask = 0;

  return;
}



///////////////////////////////////////////////////////////////////////
// create clones of the erec histograms for calculting binned FOM
void optimusPrime::deleteHistos(){

  for (int i=0; i<10; i++){
    cout<<"optimusPrime::deleteHistos: Deleting histogram "<<i<<endl;
    if(hErec[i]!=NULL) hErec[i]->Delete();
  }

  //
  return;
}



///////////////////////////////////////////////////////////////////////
// create clones of the erec histograms for calculting binned FOM
void optimusPrime::initHistos(int nselection,int oscpar){

 // get seed based on selection
 TH1D* hseed;
 if (nselection==2){
   if (oscpar==0) hseed = (TH1D*)hErecSeedNuMuTh23->Clone("hseed");
   if (oscpar==1) hseed = (TH1D*)hErecSeedNuMudM2->Clone("hseed");
 }
 else if (nselection==1){
   hseed = (TH1D*)hErecSeedNuE->Clone("hseed");
 }
 else if (nselection==3){
   hseed = (TH1D*)hErecSeedNuE->Clone("hseed");
 }
 hseed->Reset();

 // clone from seed
 hseed->SetStats(0);
 hseed->SetBit(TH1::kNoTitle,0);
 //
 hErec[0] = (TH1D*)hseed->Clone("herec_power");
 hErec[0]->SetTitle("(dN/d#theta)^{2}");
 //
 hErec[1] = (TH1D*)hseed->Clone("herec_N");
 hErec[1]->SetTitle("# of Events");
 //
 hErec[2] = (TH1D*)hseed->Clone("herec_syst");
 hErec[2]->SetTitle("Systematic Estimate");
 //
 hErec[3] = (TH1D*)hseed->Clone("herec_fom");
 hErec[3]->SetTitle("Figure of Merit");
 //
 hErec[4] = (TH1D*)hseed->Clone("herec_bg");
 hErec[4]->SetTitle("|dN/d#theta|");
 //
 hErec[5] = (TH1D*)hseed->Clone("herec_sig");
 hErec[5]->SetTitle("F.O.M.");
 //
 hErec[6] = (TH1D*)hseed->Clone("herec_ccqe");
 hErec[6]->SetTitle("CCQE");
 //
 hErec[7] = (TH1D*)hseed->Clone("herec_ccnqe");
 hErec[7]->SetTitle("CCnQE");
 //
 hErec[8] = (TH1D*)hseed->Clone("herec_ccwrong");
 hErec[8]->SetTitle("CCWrong");
 //
 hErec[9] = (TH1D*)hseed->Clone("herec_nc");
 hErec[9]->SetTitle("NC");

 // reset all bin contents
 for (int ih=0; ih<6; ih++){
    hErec[ih]->Reset();
 }
 
 hseed->Delete();

 //
 return;
}



//////////////////////////////////////////
// constructor
optimusPrime::optimusPrime(TChain* t2kmc, int nevts,
                          const char* datadir,
                          const char* nuefile,
                          const char* numufile,
                          const char* nue1rpifile){

 // set chain pointer
 chmc = t2kmc;

 // set file names of uncertainty maps
 UncMapFileNameNuE = nuefile;
 UncMapFileNameNuMu = numufile;
 UncMapFileNameNuE1RPi = nue1rpifile;

 // set event reader
 mcevent = new fqProcessedEvent(chmc);

 // max # of MC events to use
 nevents = nevts;

 // flg to toggle different variables to return
 // 0 -> figure of merit
 // 1 -> dN/dtheta
 // 2 -> (dN/dtheta) / N
 // 3 -> N
 // 4 -> sys / N
 FOMType = 0;
 
 if (nevents>chmc->GetEntries()){
   nevents = chmc->GetEntries();
   flgUseEventList = 1;
 }
 else{
  flgUseEventList = 1;
 }

 // turn off some branches for faster array filling
 chmc->SetBranchStatus("*",0);
 chmc->SetBranchStatus("mode",1);
 chmc->SetBranchStatus("wallv",1);
 chmc->SetBranchStatus("towallv",1);
 chmc->SetBranchStatus("nhitac",1);
 chmc->SetBranchStatus("attribute",1);
 chmc->SetBranchStatus("oscpower",1);
 chmc->SetBranchStatus("evtweight",1);
 chmc->SetBranchStatus("fq*",1);
 chmc->SetBranchStatus("ipnu",1);
 chmc->SetBranchStatus("pmomv",1);

 // initialze some histos
 hFVAll = new TH2FV("hall",-1,30,0,800,30,0,800);
 hFVAvg = new TH2FV("havg",-1,30,0,800,30,0,800);
 hErecSeedNuMuTh23 = new TH1D("hErecSeedNuMuTh23","hErecSeedNuMuTh23",EnuNBinsMuonTh23,EnuBinningMuonTh23);
 hErecSeedNuMudM2 = new TH1D("hErecSeedNuMudM2","hErecSeedNuMudM2",EnuNBinsMuondM2,EnuBinningMuondM2);
 hErecSeedNuE  = new TH1D("hErecSeedNuE","hErecSeedNuE",EnuNBinsElectron,EnuBinningElectron);

 // set up object to read uncertainties for each event
 uncNuE = new moreUncertainties(datadir,UncMapFileNameNuE.Data());
 uncNuE1RPi = new moreUncertainties(datadir,UncMapFileNameNuE1RPi.Data());
 uncNuMu = new moreUncertainties(datadir,UncMapFileNameNuMu.Data());

 // overall scaling factors
 Scale = 1.;
 SysScale = 1.;

 // fill large MC array for faster reading later
 fillArray();

 // use full spectrum
 flgUseSpectrum = 1;

 // init histso
 initHistos(1);

 // dont fill summary plots unless told to do so
 flgPrintSummary = 0;

 // dont use a mask unless told to do so
 flgUseMask = 0;

 // dont use floor unless told
 flgUseFloor = 0;

 // count ccnqe as signal?
 flgUseCCnQE = 0;

 // initialize summary plots (for comparisons, etc.)
 plots1 = new summaryPlots("plots_1");
 plots1->setLargeArray(fastevents);
 plots1->Init();
 //
 plots2 = new summaryPlots("plots_2");
 plots2->setLargeArray(fastevents);
 plots2->Init();

}



/////////////////////////////////////////
// get any additional uncertainty
//float optimusPrime::getMoreUncertainty(float wallv, float wallrc, float towallrc, float erec){
//  return uncertaintyCalculator->getTotalUncertainty(wallv,wallrc,towallrc,erec);
//}
/////////////////////////////////////////
// get the oscillation power for this event
float optimusPrime::getOscPower(int nutype, int oscpar){

 // NC events do not contribute 
 if (TMath::Abs(mcevent->mode)>=30) return 0.;
 
 // other nu do not contribute
 if (TMath::Abs(mcevent->ipnu[0]!=nutype)) return 0.;

 float oscpow = (float)mcevent->evtweight*(float)mcevent->oscpower[oscpar]; 

 return oscpow;

}



/////////////////////////////////////////
// get the oscillation power for this event
float optimusPrime::getOscPowerFast(int nutype, int ientry, int oscpar){


    // is dead region?
    if (fastevents->vwallv[ientry] < 0.) { return 0.;}

    // is NC?
    if (TMath::Abs(fastevents->vmode[ientry])>=30){return 0.;}

    // is Mis ID?
    if (TMath::Abs(fastevents->vnutype[ientry])!= nutype) {return 0;}

    // is CCQE?
    if (TMath::Abs(fastevents->vmode[ientry])<=10){
      float opow = getEventWeight(ientry)*fastevents->voscpower[ientry][oscpar];
      if (FOMType==2) opow = getEventWeight(ientry);
      return opow;
    }

    // CCnQE 
    if (TMath::Abs(fastevents->vmode[ientry])<30) {
      if (!flgUseCCnQE) return 0.;
      else{
        float opow = getEventWeight(ientry)*fastevents->voscpower[ientry][oscpar];
        if (FOMType==2) opow = getEventWeight(ientry);
        return opow;   
      }
    }

    // CCnQE 
    return 0.;
}



void optimusPrime::fillFVHistoFast(){
  for (int ievt=0; ievt<nevents; ievt++){
    hFVAll->Fill(fastevents->vfqtowall[ievt],fastevents->vfqwall[ievt],fastevents->vweight[ievt]);
    if (AvgType==0){
       // Signal
       if ((fastevents->vmode[ievt]==1)&&(fastevents->vnutype[ievt]==14)){
         hFVAvg->Fill(fastevents->vfqtowall[ievt],fastevents->vfqwall[ievt],fastevents->vweight[ievt]);
       }
    }
    if (AvgType==1){
       // BG 
       if (!((fastevents->vmode[ievt]==1)&&(fastevents->vnutype[ievt]==14))){
         hFVAvg->Fill(fastevents->vfqtowall[ievt],fastevents->vfqwall[ievt],fastevents->vweight[ievt]);
       } 
    }
    if (AvgType==2){    
      // Power 
      hFVAvg->Fill(fastevents->vfqtowall[ievt],fastevents->vfqwall[ievt],fastevents->vweight[ievt]*fastevents->voscpower[ievt][0]);      
    }
  }
  hFVAll->Draw("colz");
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////
void optimusPrime::calcDeltaMapNuE(float twbest, float wbest, float twmax, float wmax, int npts){

 // initialize
 double fommax = -1.;
 int    maxbin = -1;

 hDelta = new TH2FV("hdelta",-1,npts,0,twmax,npts,0,wmax);
 hDelta->SetContour(50);
 hDeltaBg = new TH2FV("hdeltabg",-1,npts,0,twmax,npts,0,wmax);
 hDeltaBg->SetContour(50);
 hDeltaSg = new TH2FV("hdeltasg",-1,npts,0,twmax,npts,0,wmax);
 hDeltaSg->SetContour(50);

 int nbincount = 0;
 for (int ibin = 0; ibin<hDelta->GetNumberOfBins(); ibin++){

   float wallcut = (float)hDelta->GetBinCenterY(ibin);
   float towallcut = (float)hDelta->GetBinCenterX(ibin);

   if (towallcut<=wallcut) continue;
   float value = calcDeltaNuE(twbest,wbest,towallcut,wallcut);
   hDelta->SetBinContent(ibin,TMath::Abs(DeltaSg)+TMath::Abs(DeltaBg));
   hDeltaSg->SetBinContent(ibin,DeltaSg);
   hDeltaBg->SetBinContent(ibin,DeltaBg);

 }
 hDelta->SetMinimum(0);
 hDelta->SetMaximum(1.);
 hDelta->Draw("colz");
 return;
  
}




///////////////////////////////////////////////////////////////////////////////////////////////////////////
void optimusPrime::calcDeltaMapNuMu(float twbest, float wbest, float twmax, float wmax, int npts){

 // initialize
 double fommax = -1.;
 int    maxbin = -1;

 hDelta = new TH2FV("hdelta",-1,npts,0,twmax,npts,0,wmax);
 hDelta->SetContour(50);
 hDeltaBg = new TH2FV("hdeltabg",-1,npts,0,twmax,npts,0,wmax);
 hDeltaBg->SetContour(50);
 hDeltaSg = new TH2FV("hdeltasg",-1,npts,0,twmax,npts,0,wmax);
 hDeltaSg->SetContour(50);

 int nbincount = 0;
 for (int ibin = 0; ibin<hDelta->GetNumberOfBins(); ibin++){

   float wallcut = (float)hDelta->GetBinCenterY(ibin);
   float towallcut = (float)hDelta->GetBinCenterX(ibin);

   if (towallcut<=wallcut) continue;
   float value = calcDeltaNuMu(twbest,wbest,towallcut,wallcut);

   hDelta->SetBinContent(ibin,TMath::Abs(DeltaSg)+TMath::Abs(DeltaBg));
   hDeltaSg->SetBinContent(ibin,DeltaSg);
   hDeltaBg->SetBinContent(ibin,DeltaBg);

 }
 hDelta->SetMinimum(0);
 hDelta->SetMaximum(1.);
 hDelta->Draw("colz");
  return;
  
}


//////////////////////////////////////////////////////////////////////////////////////////////////
// useful method to map out the figure of merit in each bin
void optimusPrime::calcFOMMap(float towallmax, float wallmax,int oscpar, int npts, int nselection){

 // initialize
 double fommax = -1.;
 int    maxbin = -1;

 // initialize FV bin map
 hFV = new TH2FV("h",-1,npts,0,towallmax,npts,0,wallmax);
 hFV->SetContour(50);

 // delete histograms
 deleteHistos();

 // create erec histos depending on the sample you want
 initHistos(nselection,oscpar);

 // loop over FV bins, for each bin calculate the figure of merit
 for (int ibin = 0; ibin<hFV->GetNumberOfBins(); ibin++){

   // get the wall cuts from the center of the current bin
   float wallcut = (float)hFV->GetBinCenterY(ibin);
   float towallcut = (float)hFV->GetBinCenterX(ibin);

   // reject the impossible
   if (towallcut<=wallcut) continue;

   // apply mask?
   // can use another TH2FV to select only a few bins for FOM evaluation
   if (flgUseMask){
     int maskbin = hMask->FindBin(towallcut,wallcut);
     double maskval = hMask->GetBinContent(maskbin);
     if (maskval<maskThresh) continue;
   }

   // find the FOM value for this set of cuts
   float value = calcFOMBinned(nselection, towallcut, wallcut, oscpar, 1);

   // keep track of the best point
   if (value>fommax){
     fommax = value;
     maxbin = ibin;
   }

   // fill the map
   hFV->SetBinContent(ibin,value);
 }

 // print some results
 cout<<"Max value: "<<fommax<<endl;
 cout<<"Max bin: "<<maxbin<<endl;
 bestFOMbin = maxbin;
 bestFOMvalue = fommax;
 cout<<"Best Towall: "<<hFV->GetBinCenterX(maxbin)<<endl;
 cout<<"Best Wall: "<<hFV->GetBinCenterY(maxbin)<<endl;
 hFV->SetMinimum(0.7*fommax);
 hFV->SetMaximum(1.01*fommax);
 
 // draw result
// hFV->Draw("colz");
 drawFOM(nselection, oscpar);
 //
 return;

}



//////////////////////////////////////////
//Read events into memory for fast looping
void optimusPrime::fillArray(){

  fastevents = new mcLargeArray(chmc,nevents); 

  //
  return;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////
float optimusPrime::calcFOMBinned(int nselection, float towallmin, float wallmin, int oscpar, int iplt){

  // summarize?
  if (flgPrintSummary){
    if (iplt==1){
      plots1->clearHistos();
    }
    if (iplt==2){
      plots2->clearHistos();
    }
  }

  // towall must be smaller than wall
  if (towallmin<wallmin){
    return 0.;
  }

  // reject if below floor
  if (flgUseFloor){
    if (towallmin<150.){
      return 0.;
    }
    if (wallmin<50.){
      return 0.;
    }
  }

  // reset values
  Nevents = 0.;
  NB = 0.; 
  NS = 0.; 
  Power = 0.;
  Syst = 0.;
 
  // clear previous histos
  hErec[0]->Reset();
  hErec[1]->Reset();
  hErec[2]->Reset();
  hErec[3]->Reset();
  hErec[4]->Reset();

  // make arrays
  const int nn = hErec[0]->GetNbinsX()+1;
  float Pow[nn];
  float Nev[nn];
  float Sys[nn];
  for (int i=0; i<nn; i++){
    Pow[i]=0.;
    Nev[i]=0.;
    Sys[i]=0.;
  }

  // neutrino type of the selection you desire
  int nu_type = 0;
  if (nselection==1){
    nu_type = 12;
    uncertaintyCalculator = uncNuE;
    flgUseCCnQE = 0;
  }
  else if (nselection==2){
    nu_type = 14;
    uncertaintyCalculator = uncNuMu;
    flgUseCCnQE = 0;
  }
  else if (nselection==3){
    nu_type = 12;
    uncertaintyCalculator = uncNuE1RPi;
    flgUseCCnQE = 1;
  }

  // loop over events
  for (int i=0; i<nevents; i++){

     // see if event passes your selection cuts
     int ipass = 0.;
     if (nselection==1){
       ipass = fastevents->vpassnue[i];
     }
     else if (nselection==2){
       ipass = fastevents->vpassnumu[i];
     }
     else if (nselection==3){
       ipass = fastevents->vpassnue1rpi[i];
     }

     // if it passes, add to spectrum
     if (ipass){

      // ..as long as it passes the FV cuts
      if ((fastevents->vfqwall[i] >= wallmin)&&(fastevents->vfqtowall[i]>=towallmin)){
       
        // get erec bin of this event and 
        int erecbin = 0.;
        if (nselection==1 || nselection==3){
          erecbin = hErec[0]->FindBin(fastevents->vfqenue[i]);
        }
        else if (nselection==2){
          erecbin = hErec[0]->FindBin(fastevents->vfqenumu[i]);
        }
       
        // get oscillation power (derivative w.r.t. oscillation parameter)
        float opow = getOscPowerFast(nu_type,i,oscpar)*Scale;
        Pow[erecbin] += opow;
       
        // keep track of all events
        float ww = getEventWeight(i)*Scale;
        Nev[erecbin] += ww;
       
        // keep track of systematic uncertainty
        float syserr =  getSystUncertainty(i,nu_type)*Scale*SysScale;
        Sys[erecbin] += syserr;

        // plots for the details
        if (flgPrintSummary){
          if (iplt==1) plots1->fillAllFromArray(i,opow,syserr);
          if (iplt==2) plots2->fillAllFromArray(i,opow,syserr);
        }

      }
    }
  }

  // add up figure of merit in each bin
  return calcFOM(Pow,Nev,Sys,nn);


}



/////////////////////////////////////////////
// draw stacked histogram of events
void optimusPrime::showBreakdown(){

  hs = new THStack("hs","");
  
  // color scheme
  int colerz[5] = {4,7,2,1,6};

  hErec[6]->SetFillColor(kCyan);
  hErec[7]->SetFillColor(kBlue);
  hErec[8]->SetFillColor(kRed);
  hErec[9]->SetFillColor(kBlack);

  hs->Add(hErec[9]);
  hs->Add(hErec[8]);
  hs->Add(hErec[7]);
  hs->Add(hErec[6]);

  hs->Draw("h");
  
  return;
}



///////////////////////////////////////////////////////
// get the weight for this event (can increase norm of
// specifice event catagories i.e. entering b.g.
float optimusPrime::getEventWeight(int iev){

  float ww = 1.0;

  ww *= fastevents->vweight[iev];

  if (fastevents->vwallv[iev] < 0.){
    ww *= 1.1; //< 10% additional weight 
  }

  return ww;
}



/////////////////////////////////////////////////////
double optimusPrime::calcDeltaNuE(float tw1, float w1, float tw2, float w2){

  // clear histograms
  plots1->clearHistos();
  plots2->clearHistos();

  // loop over events
  for (int i=0; i<nevents; i++){

     // see if event passes numu cuts
     int ipass = fastevents->vpassnue[i];

     // if it passes, add to spectrum
     if (ipass){

      // what about FV cuts?
      int passfvcuts1 = 0;
      int passfvcuts2 = 0;
      if ((fastevents->vfqwall[i] >= w1)&&(fastevents->vfqtowall[i]>=tw1)) passfvcuts1 = 1;
      if ((fastevents->vfqwall[i] >= w2)&&(fastevents->vfqtowall[i]>=tw2)) passfvcuts2 = 1;

      // if passes, move on
      if ( passfvcuts1 || passfvcuts2){
        
        // fill histos
        if (passfvcuts1) plots1->fillAllFromArray(i,0.,0.);
        if (passfvcuts2) plots2->fillAllFromArray(i,0.,0.);

      }
    }
  }

  // print absolute diff
  Delta = getAbsDifference(plots1->pltEnuE,plots2->pltEnuE);
  DeltaSg = getAbsDifference(plots1->pltEnuESg,plots2->pltEnuESg);
  DeltaBg = getAbsDifference(plots1->pltEnuEBg,plots2->pltEnuEBg);

  return Delta;

}



/////////////////////////////////////////////////////
double optimusPrime::calcDeltaNuMu(float tw1, float w1, float tw2, float w2){

  // clear histograms
  plots1->clearHistos();
  plots2->clearHistos();

  // loop over events
  for (int i=0; i<nevents; i++){

     // see if event passes numu cuts
     int ipass = fastevents->vpassnumu[i];

     // if it passes, add to spectrum
     if (ipass){

      // what about FV cuts?
      int passfvcuts1 = 0;
      int passfvcuts2 = 0;
      if ((fastevents->vfqwall[i] >= w1)&&(fastevents->vfqtowall[i]>=tw1)) passfvcuts1 = 1;
      if ((fastevents->vfqwall[i] >= w2)&&(fastevents->vfqtowall[i]>=tw2)) passfvcuts2 = 1;

      // if passes, move on
      if ( passfvcuts1 || passfvcuts2){
        
        // fill histos
        if (passfvcuts1) plots1->fillAllFromArray(i,0.,0.);
        if (passfvcuts2) plots2->fillAllFromArray(i,0.,0.);

      }
    }
  }

  // print absolute diff
  Delta = getAbsDifference(plots1->pltEnuMu,plots2->pltEnuMu);
  DeltaSg = getAbsDifference(plots1->pltEnuMuSg,plots2->pltEnuMuSg);
  DeltaBg = getAbsDifference(plots1->pltEnuMuBg,plots2->pltEnuMuBg);

  return Delta;

}


/////////////////////////////////////////////////////
int optimusPrime::isSmallDifference(float tw1, float w1, float tw2, float w2, int oscpar, int nselection){

  // make sure to record histograms
  flgPrintSummary = 1;

  // do the calculations and fill plots1
  float fom1 = 0.;
  float fom2 = 0.;

  fom1 = calcFOMBinned(nselection,tw1,w1,oscpar,1);
  fom1 = calcFOMBinned(nselection,tw2,w2,oscpar,2);

  double chi2=0.;
  double ksP=0.;
  double absdiff=0.;

  if (nselection==2){
    plots1->pltEnuMu->SetLineColor(kBlue);
    plots2->pltEnuMu->SetLineColor(kRed);
    plots1->pltEnuMu->Draw("e");
    plots2->pltEnuMu->Draw("same");
    chi2 = plots1->pltEnuMu->Chi2Test(plots2->pltEnuMu,"WW");
    ksP = plots1->pltEnuMu->KolmogorovTest(plots2->pltEnuMu,"N");
    absdiff = getAbsDifference(plots1->pltEnuMu,plots2->pltEnuMu );
  }
  else if (nselection==1 || nselection==3){
    plots1->pltEnuE->SetLineColor(kBlue);
    plots2->pltEnuE->SetLineColor(kRed);
    plots1->pltEnuE->Draw("e");
    plots2->pltEnuE->Draw("same");
    chi2 = plots1->pltEnuE->Chi2Test(plots2->pltEnuE,"WW");
    ksP = plots1->pltEnuE->KolmogorovTest(plots2->pltEnuE,"N");
  }

  return 1;
}



/////////////////////////////////////////////////////
// compare tow sets of FV cuts
void optimusPrime::compareFOM(float tw1, float w1, float tw2, float w2, int oscpar, int nselection){
  
  // make sure to record histograms
  flgPrintSummary = 1;

  // canvas setup
  multiPad = new TCanvas("multiPad","multiPad",700,800);
  multiPad->Divide(2,2);

  // colerz
  int colerz[5] = {4,7,2,15,6};

  // set erec binning
  deleteHistos();
  initHistos(nselection,oscpar);

  // do the calculations and fill plots1
  float fom1 = 0.;
  fom1 = calcFOMBinned(nselection,tw1,w1,oscpar,1);

  TH1D* htmp[4];
  htmp[0]  = (TH1D*)hErec[0]->Clone("htmp1");
  htmp[1]  = (TH1D*)hErec[1]->Clone("htmp2");
  htmp[2]  = (TH1D*)hErec[2]->Clone("htmp3");
  htmp[3]  = (TH1D*)hErec[3]->Clone("htmp4");

  // new cut (1)
  multiPad->cd(1);
  htmp[0]->SetLineStyle(1);
  htmp[0]->SetLineWidth(3);
  htmp[0]->SetLineColor(1);
  htmp[0]->Draw();
  multiPad->cd(2);
  htmp[1]->SetLineStyle(1);
  htmp[1]->SetLineWidth(3);
  htmp[1]->SetLineColor(1);
  htmp[1]->Draw();
  multiPad->cd(3);
  htmp[2]->SetLineStyle(1);
  htmp[2]->SetLineWidth(3);
  htmp[2]->SetLineColor(1);
  htmp[2]->SetMinimum(0);
  htmp[2]->SetMaximum( htmp[1]->GetMaximum() );
  htmp[2]->Draw();
  multiPad->cd(4);
  htmp[3]->SetLineStyle(1);
  htmp[3]->SetLineWidth(3);
  htmp[3]->SetLineColor(1);
  htmp[3]->Draw();

  // do the calculations and fill plots2
  float fom2 = 0.;
  fom1 = calcFOMBinned(nselection,tw2,w2,oscpar,2);

  // old cut (2)
  multiPad->cd(1);
  hErec[0]->SetLineStyle(2);
  hErec[0]->SetLineWidth(3);
  hErec[0]->SetLineColor(2);
  hErec[0]->Draw("same");
  multiPad->cd(2);
  hErec[1]->SetLineStyle(2);
  hErec[1]->SetLineWidth(3);
  hErec[1]->SetLineColor(2);
  hErec[1]->Draw("same");
  multiPad->cd(3);
  hErec[2]->SetLineStyle(2);
  hErec[2]->SetLineWidth(3);
  hErec[2]->SetLineColor(2);
  hErec[2]->Draw("same");
  multiPad->cd(4);
  hErec[3]->SetLineStyle(2);
  hErec[3]->SetLineWidth(3);
  hErec[3]->SetLineColor(2);
  hErec[3]->Draw("same");
  
  //
  cout<<"FOM1: "<<fom1<<endl;
  cout<<"FOM2: "<<fom2<<endl;

  cout<<"KS N: "<<htmp[1]->KolmogorovTest(hErec[1],"N")<<endl;
  cout<<"KS S: "<<htmp[2]->KolmogorovTest(hErec[2],"N")<<endl;

  // turn recording back off
  flgPrintSummary = 0;

  //
  return;
}


/////////////////////////////////////////////////////
// compare tow sets of FV cuts
void optimusPrime::compareCuts(float tw1, float w1, float tw2, float w2, int oscpar, int nselection){
  
  // make sure to record histograms
  flgPrintSummary = 1;

  // canvas setup
  multiPad = new TCanvas("multiPad","multiPad",700,800);
  multiPad->Divide(2,3);

  // set erec binning
  deleteHistos();
  initHistos(nselection,oscpar);

  // colerz
  int colerz[5] = {4,7,2,15,6};

  // do the calculations and fill plots1
  float fom1 = 0.;

  fom1 = calcFOMBinned(nselection,tw1,w1,oscpar,1);

  // do the calculations and fill plots2
  float fom2= calcFOMBinned(nselection,tw2,w2,oscpar,2);


  if (nselection==2){ 
    multiPad->cd(1);
    plots1->pltEnuMu->Draw();
    plots1->pltEnuMu->SetLineWidth(3);
    plots2->pltEnuMu->SetLineWidth(3);
    plots2->pltEnuMu->SetLineStyle(2);
    plots1->pltEnuMu->Draw("h");
    plots2->pltEnuMu->Draw("sameh");
    double ymin = 0.;
    double ymax = plots1->pltEnuMu->GetMaximum();
    for (int icat=0; icat<5; icat++){
      plots1->pltEnuMuCat[icat]->SetLineColor(colerz[icat]);
      plots2->pltEnuMuCat[icat]->SetLineColor(colerz[icat]);
      plots2->pltEnuMuCat[icat]->SetLineWidth(3);
      plots1->pltEnuMuCat[icat]->SetLineWidth(3);
      plots2->pltEnuMuCat[icat]->SetLineStyle(2);
      multiPad->cd(icat+2);
      plots1->pltEnuMuCat[icat]->Draw("h");
      plots2->pltEnuMuCat[icat]->Draw("sameh");

    }
  }
  else if (nselection==1){
    multiPad->cd(1);
    plots1->pltEnuE->Draw();
    plots1->pltEnuE->SetLineWidth(3);
    plots2->pltEnuE->SetLineWidth(3);
    plots2->pltEnuE->SetLineStyle(2);
    plots1->pltEnuE->Draw("h");
    plots2->pltEnuE->Draw("sameh");
    double ymin = 0.;
    double ymax = plots1->pltEnuE->GetMaximum();
    for (int icat=0; icat<5; icat++){
      plots1->pltEnuECat[icat]->SetLineColor(colerz[icat]);
      plots2->pltEnuECat[icat]->SetLineColor(colerz[icat]);
      plots2->pltEnuECat[icat]->SetLineWidth(3);
      plots1->pltEnuECat[icat]->SetLineWidth(3);
      plots2->pltEnuECat[icat]->SetLineStyle(2);
      multiPad->cd(icat+2);
      plots1->pltEnuECat[icat]->Draw("h");
      plots2->pltEnuECat[icat]->Draw("sameh");
    }
  }
  else if (nselection==3){
    multiPad->cd(1);
    plots1->pltEnuE1RPi->Draw();
    plots1->pltEnuE1RPi->SetLineWidth(3);
    plots2->pltEnuE1RPi->SetLineWidth(3);
    plots2->pltEnuE1RPi->SetLineStyle(2);
    plots1->pltEnuE1RPi->Draw("h");
    plots2->pltEnuE1RPi->Draw("sameh");
    double ymin = 0.;
    double ymax = plots1->pltEnuE1RPi->GetMaximum();
    for (int icat=0; icat<5; icat++){
      plots1->pltEnuE1RPiCat[icat]->SetLineColor(colerz[icat]);
      plots2->pltEnuE1RPiCat[icat]->SetLineColor(colerz[icat]);
      plots2->pltEnuE1RPiCat[icat]->SetLineWidth(3);
      plots1->pltEnuE1RPiCat[icat]->SetLineWidth(3);
      plots2->pltEnuE1RPiCat[icat]->SetLineStyle(2);
      multiPad->cd(icat+2);
      plots1->pltEnuE1RPiCat[icat]->Draw("h");
      plots2->pltEnuE1RPiCat[icat]->Draw("sameh");
    }
  }

  //
  cout<<"FOM1: "<<fom1<<endl;
  cout<<"FOM2: "<<fom2<<endl;

  // turn recording back off
  flgPrintSummary = 0;

  //
  return;
}



///////////////////////////////////////////////////
// calculate FOM from arrays 
float optimusPrime::calcFOM(float* pow, float* nev, float* sys, int nbin){

 if (flgPrintSummary){
   hErec[0]->Reset();
   hErec[1]->Reset();
   hErec[2]->Reset();
   hErec[3]->Reset();
   hErec[4]->Reset();
 }

 float fom_total = 0.; 
 float syst_total  = 0.;
 float power_total = 0.;
 float nev_total = 0.;

 // loop over spectrum bins
 for (int i=0; i<nbin; i++){
  
    // figure of merit in this spectrum bin
    float fom_thisbin = 0.;

    // if there are events in this bin, it contributes to total
    if (nev[i]>0.){
      // optimization figure of merit
      if (FOMType==0) fom_thisbin = (pow[i]*pow[i])/((sys[i]*sys[i])+nev[i]);
      // other pototentially interesting things
      if (FOMType==1) fom_thisbin = (pow[i]);
      if (FOMType==2){
        if (nev[i]!=0.) fom_thisbin = (pow[i])/(nev[i]);
      }
      if (FOMType==3) fom_thisbin = (nev[i]);
      if (FOMType==4){
        if (nev[i]!=0.) fom_thisbin = (sys[i]/nev[i]);
      }
    }
   
    if (flgPrintSummary){
      hErec[0]->SetBinContent(i,(pow[i]*pow[i]));
      hErec[1]->SetBinContent(i,nev[i]);
      hErec[1]->SetBinError(i,TMath::Sqrt(nev[i]));
      hErec[2]->SetBinContent(i,sys[i]);
      hErec[3]->SetBinContent(i,fom_thisbin);
      hErec[4]->SetBinContent(i,TMath::Abs(pow[i]));
    }

    fom_total += fom_thisbin;
    power_total += (pow[i]*pow[i]);
    syst_total  += (sys[i]);
    nev_total   += nev[i];
 }

 // print some info?
 if (flgPrintSummary){
   cout<<"total P: "<<power_total<<endl;
   cout<<"total N: "<<nev_total<<endl;
   cout<<"total S: "<<syst_total<<endl;
   cout<<"FOM: "<<fom_total<<endl;
   smallVariation = defineSmall2(nbin,pow,nev,sys);
   cout<<"Small: "<<smallVariation<<endl;;
 }

// if (FOMType==1) return power_total;
// if (FOMType==2) return nev_total;
// if (FOMType==3) return syst_total;
// if (FOMType==4) return power_total/nev_total;
// if (FOMType==5) return syst_total/nev_total;

// if (FOMType==3) return syst_total/TMath::Sqrt(nev_total);
// if (FOMType==4) return fom_total/(syst_total*syst_total);
// if (FOMType==5) return fom_total/(syst_total);
// if (FOMType==6) return fom_total;
 if (FOMType==8) return fom_total;
 return fom_total;

}


//////////////////////////////////////////
// calculate some interesting FV maps
//  1) number of events in each FV bin
//  2) power in each FV bin
//  3) systematic error in each FV bin
//  4) fom in each FV bin
void optimusPrime:: calcFVSummary(int oscpar, int nutype){

 // make histograms
 int nbins = 30;
 double towallmax = 1200;
 double wallmax = 1200;
 hFVSummary[0] = new TH2FV("hfv_nev",-1,nbins,0,towallmax,nbins,0,wallmax);
 hFVSummary[1] = new TH2FV("hfv_pow",-1,nbins,0,towallmax,nbins,0,wallmax);
 hFVSummary[2] = new TH2FV("hfv_syst",-1,nbins,0,towallmax,nbins,0,wallmax);
 hFVSummary[3] = new TH2FV("hfv_fom",-1,nbins,0,towallmax,nbins,0,wallmax);
 hFVSummary[4] = new TH2FV("hfv_enutrue",-1,nbins,0,towallmax,nbins,0,wallmax);
 hFVSummary[5] = new TH2FV("hfv_enurc",-1,nbins,0,towallmax,nbins,0,wallmax);
 hFVSummary[6] = new TH2FV("hfv_ccqe",-1,nbins,0,towallmax,nbins,0,wallmax);
 hFVSummary[7] = new TH2FV("hfv_nccqe",-1,nbins,0,towallmax,nbins,0,wallmax);
 // 1D
 hSummary[0] = new TH1D("hmode","hmode",160,-40,40);
 hSummary[1] = new TH1D("henuv","henuv",100,0,5000);
 hSummary[2] = new TH1D("hwall","hwall",100,0,3000);
 hSummary[3] = new TH1D("ht0wall","htowall",100,0,3000);
 hSummary[4] = new TH1D("hwallv","hwallv",100,-500,3000);
 hSummary[5] = new TH1D("nring","nring",10,0,10);
 hSummary[6] = new TH1D("pid","pid",100,-3000,3000);


 // fill histograms
 for (int iev=0; iev<nevents; iev++){
   if (nutype==14){
     if (!passNuMuCuts(iev)) continue;
   }
   else{
     if (!passNuECuts(iev)) continue;
   }
   // fill 
   float ww = fastevents->vweight[iev];
   float sys = getSystUncertainty(iev,nutype); 
   float pow = getOscPowerFast(nutype,iev,oscpar);
   float enurc = 0;
   if (nutype==14) enurc = fastevents->vfqenumu[iev];
   else { enurc = fastevents->vfqenumu[iev];}
   float enuv = fastevents->vpmomv[iev];
   // use RC
   hFVSummary[0]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],ww);
   hFVSummary[1]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],pow);
   hFVSummary[2]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],sys);
   hFVSummary[4]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],enuv*ww);
   hFVSummary[5]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],enurc*ww);
   // use true
//   hFVSummary[0]->Fill(fastevents->vtowallv[iev],fastevents->vwallv[iev],ww);
//   hFVSummary[1]->Fill(fastevents->vtowallv[iev],fastevents->vwallv[iev],pow);
//   hFVSummary[2]->Fill(fastevents->vtowallv[iev],fastevents->vwallv[iev],sys);
//   hFVSummary[4]->Fill(fastevents->vtowallv[iev],fastevents->vwallv[iev],enuv*ww);
//   hFVSummary[5]->Fill(fastevents->vtowallv[iev],fastevents->vwallv[iev],enurc*ww);

   if (TMath::Abs(fastevents->vmode[iev])<30) hFVSummary[6]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],ww);
   else {hFVSummary[7]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],ww);}

   hSummary[0]->Fill(fastevents->vmode[iev],ww);
   hSummary[1]->Fill(fastevents->vpmomv[iev],ww);
   hSummary[2]->Fill(fastevents->vfqwall[iev],ww);
   hSummary[3]->Fill(fastevents->vfqtowall[iev],ww);
   hSummary[4]->Fill(fastevents->vwallv[iev],ww);
   hSummary[5]->Fill(fastevents->vfqnring[iev],ww);
   hSummary[6]->Fill(fastevents->vfqpid[iev],ww);
 }

 // calculate fom histogram and normalize
 for (int ibin=1; ibin<hFVSummary[0]->GetNumberOfBins(); ibin++){
   float nev = hFVSummary[0]->GetBinContent(ibin);
   if (nev==0) continue;
   float pow = hFVSummary[1]->GetBinContent(ibin);
   float sys = hFVSummary[2]->GetBinContent(ibin);
   float fom = (pow*pow)/((sys*sys) + nev);
   hFVSummary[1]->SetBinContent(ibin,TMath::Abs(pow)/nev);
   hFVSummary[2]->SetBinContent(ibin,sys/nev);
   hFVSummary[3]->SetBinContent(ibin,fom/nev);
   hFVSummary[4]->SetBinContent(ibin, hFVSummary[4]->GetBinContent(ibin)/nev);
   hFVSummary[5]->SetBinContent(ibin, hFVSummary[5]->GetBinContent(ibin)/nev);
 }

 hFVSummary[3]->Draw("colz");

 return;
}



////////////////////////////////////////////////////////////
// calculate the fractional weighted systematic uncertainty
// for an event
float optimusPrime::getSystUncertainty(int i, int nutype){

   // is the lepton mis-IDed?
   int wronglepton = 0;
   if (TMath::Abs(fastevents->vnutype[i])!=nutype) wronglepton = 1;

   float sys = uncertaintyCalculator->getTotalUncertainty(fastevents->vwallv[i],
                                                          fastevents->vfqwall[i],
                                                          fastevents->vfqtowall[i],
                                                          fastevents->vfqenumu[i],
                                                          fastevents->vmode[i],
                                                          wronglepton);

   return sys*getEventWeight(i);                                                       
}



/////////////////////////////////////////////////////////////////
// apply the cuts to a modified event and see if it passes
int optimusPrime::applyCutsToModifiedEvent(int iev){

  // fill tmp array with "nominal" MC values
  const int natt = 4;
  float attributesTmp[natt];
  for (int iatt=0; iatt<natt; iatt++){
    attributesTmp[iatt] = fastevents->vattribute[iev][iatt];   
  }
 
  // modify tmp array by applying the histogram shape parameters
  modifier->applyPars(fastevents->vbin[iev],
                      fastevents->vcomponent[iev],
                      attributesTmp,
                      natt);

  // fill cut parameter structure using modified attributes
  if (indexPIDPar>=0) cutPars.fqpid = attributesTmp[indexPIDPar];
  if (indexPi0Par>=0) cutPars.fqpi0par = attributesTmp[indexPi0Par];
  if (indexPiPPar>=0) cutPars.fqpippar = attributesTmp[indexPiPPar];
//  if (indexRCPar>=0) cutPars.fqrcpar = attributesTmp[indexRCPar];
  cutPars.fqrcpar = fastevents->vfqrcpar[iev];

  // other cut pars that are not modified
  cutPars.fqmome = fastevents->vfqmumom[iev];
  cutPars.fqmommu = fastevents->vfqemom[iev];
  cutPars.nhitac = fastevents->vnhitac[iev];
  cutPars.fqnsubev = fastevents->vfqnsubev[iev];
  cutPars.fqenue = fastevents->vfqenue[iev];
  cutPars.fqenumu = fastevents->vfqenumu[iev];

  // see if it passes cuts
  int passnue = selectNuE(cutPars);
  int passnumu = selectNuMu(cutPars);
  
  //
  if (passnue>0) return 1;
  if (passnumu>0) return 2;
  return 0;
  
}



///////////////////////////////////////////
// does event pass cuts?
int optimusPrime::passNuMuCuts(int i){

  int ipass = fastevents->vpassnumu[i];
  return ipass;

}



/////////////////////////////////////////
// what about nu-e cuts?
int optimusPrime::passNuECuts(int i){

                         
  int ipass = fastevents->vpassnue[i];
  return ipass;

}

#endif







