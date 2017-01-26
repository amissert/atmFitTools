#ifndef PRIME_CXX
#define PRIME_CXX


#include "optimusPrime.h"

//////////////////////////////////////////
// constructor
optimusPrime::optimusPrime(TChain* t2kmc, int nevts, const char* datadir, const char* mapfile){

 // set chain pointer
 chmc = t2kmc;

 // set event reader
 mcevent = new fqProcessedEvent(chmc);

 // max # of MC events to use
 nevents = nevts;

 // flg to toggle using energy spectrom when evaluating f.o.m
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

 // set up object to read uncertainties for each event
 uncertaintyCalculator = new moreUncertainties(datadir, mapfile);

 // overall scaling factors
 Scale = 1.;
 SysScale = 1.;

 // fill large MC array for faster reading later
 fillArray();

 // setup recon energy histos
 // seed from histogram binning in uncertainty map file
 TH1D* hseed = uncertaintyCalculator->hERecUnc[0];
// TH1D* hseed = new TH1D("hseed","hseed",1,0,8000);
 hseed->SetStats(0);
 hseed->SetBit(TH1::kNoTitle,0);
 hErec[0] = (TH1D*)hseed->Clone("herec_power");
 hErec[0]->SetTitle("Oscillation Power");
 hErec[1] = (TH1D*)hseed->Clone("herec_N");
 hErec[1]->SetTitle("# of Events");
 hErec[2] = (TH1D*)hseed->Clone("herec_syst");
 hErec[2]->SetTitle("Systematic Uncertainty");
 hErec[3] = (TH1D*)hseed->Clone("herec_fom");
 hErec[3]->SetTitle("Figure of Merit");
 hErec[4] = (TH1D*)hseed->Clone("herec_bg");
 hErec[4]->SetTitle("# BG");
 hErec[5] = (TH1D*)hseed->Clone("herec_sig");
 hErec[5]->SetTitle("F.O.M.");
 hErec[6] = (TH1D*)hseed->Clone("herec_ccqe");
 hErec[6]->SetTitle("CCQE");
 hErec[7] = (TH1D*)hseed->Clone("herec_ccnqe");
 hErec[7]->SetTitle("CCnQE");
 hErec[8] = (TH1D*)hseed->Clone("herec_ccwrong");
 hErec[8]->SetTitle("CCWrong");
 hErec[9] = (TH1D*)hseed->Clone("herec_nc");
 hErec[9]->SetTitle("NC");
 // reset all bin contents
 for (int ih=0; ih<6; ih++){
    hErec[ih]->Reset();
 }

 // use full spectrum
 flgUseSpectrum = 1;

 // dont fill summary plots unless told to do so
 flgPrintSummary = 0;

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

/*
 // NC events do not contribute 
 // Non-CCQE do not contribute
 if (TMath::Abs(fastevents->vmode[ientry])>=10) return 0.;
 
 // other nu do not contribute
 if (TMath::Abs(fastevents->vnutype[ientry]) != nutype) return 0.;

 // outside FV do not contribute
 if (fastevents->vwallv[ientry] <=0.) return 0.;

 // return power multiplied by weight
 float oscpow = fastevents->vweight[ientry]*fastevents->voscpower[ientry][oscpar]; 

 if (FOMType==2) oscpow = fastevents->vweight[ientry];

*/

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
    if (TMath::Abs(fastevents->vmode[ientry])<30) {return 0.;}

    return 0.;
}

void optimusPrime::fillFVHistoFast(){
//  hFVAll = new TH2FV("hall",-1,40,0,1000,40,0,1000);
//  hFVAvg = new TH2FV("havg",-1,40,0,1000,40,0,1000);
//  hFVAll->Reset();
//  hFVAvg->Reset();
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


//////////////////////////////////////////////////////////////////////////////////////////////////
// useful method to map out the figure of merit in each bin
void optimusPrime::calcFOMMap(float towallmax, float wallmax,int oscpar, int npts){

 // initialize
 double fommax = -1.;
 int    maxbin = -1;

 // for finding FV 
 hFV = new TH2FV("h",-1,npts,0,towallmax,npts,0,wallmax);
 hFV->SetContour(50);
 for (int ibin = 0; ibin<hFV->GetNumberOfBins(); ibin++){
   float wallcut = (float)hFV->GetBinCenterY(ibin);
   float towallcut = (float)hFV->GetBinCenterX(ibin);
   if (towallcut<=wallcut) continue;
   float value = 0.;
   if (!flgUseSpectrum) value = calcNuMuFOM(towallcut,wallcut,oscpar);
   else{ value = calcFOMSpectrumNuMu(towallcut,wallcut,oscpar);}
   if (value>fommax){
     fommax = value;
     maxbin = ibin;
   }
   hFV->SetBinContent(ibin,value);
 }
 
 hFV->Draw("colz");
 cout<<"Max value: "<<fommax<<endl;
 cout<<"Max bin: "<<maxbin<<endl;
  return;
}


//////////////////////////////////////////////////////////////////////////////////////////////////
// useful method to map out the figure of merit in each bin
void optimusPrime::calcFOMMapE(float towallmax, float wallmax,int oscpar, int npts){

 double fommax = -1.;
 int    maxbin = -1;

 hFV = new TH2FV("h",-1,npts,0,towallmax,npts,0,wallmax);
 hFV->SetContour(50);
 for (int ibin = 0; ibin<hFV->GetNumberOfBins(); ibin++){
   float wallcut = (float)hFV->GetBinCenterY(ibin);
   float towallcut = (float)hFV->GetBinCenterX(ibin);
   if (towallcut<=wallcut) continue;
   float value = calcFOMSpectrumNuE(towallcut,wallcut,oscpar);
   if (value>fommax){
     fommax = value;
     maxbin = ibin;
   }
   hFV->SetBinContent(ibin,value);
 }
 
 hFV->Draw("colz");
 cout<<"Max value: "<<fommax<<endl;
 cout<<"Max bin: "<<maxbin<<endl;
  return;
}


//////////////////////////////////////////
//Read events into memory for fast looping
void optimusPrime::fillArray(){

  fastevents = new mcLargeArray(chmc,nevents); 

  //
  return;
}


//////////////////////////////////////////////////
// get FOM using a spectrum for electron selection
float optimusPrime::calcFOMSpectrumNuE(float towallmin, float wallmin, int oscpar, int iplt){

  
  /* the old way

  // towall must be smaller than wall
  if (towallmin<wallmin){
    return 0.;
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
  hErec[5]->Reset();
  hErec[6]->Reset();
  hErec[7]->Reset();
  hErec[8]->Reset();
  hErec[9]->Reset();

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

  // loop over events
  for (int i=0; i<nevents; i++){

     // see if event passes numu cuts
     int ipass = passNuECuts(i);

     // if it passes, add to spectrum
     if (ipass){

      // ..as long as it passes the FV cuts
      if ((fastevents->vfqwall[i] >= wallmin)&&(fastevents->vfqtowall[i]>=towallmin)){
       
        // get erec bin of this event
        int erecbin = hErec[0]->FindBin(fastevents->vfqenue[i]);
        
        Pow[erecbin] += getOscPowerFast(12,i,oscpar)*Scale;
        
        Nev[erecbin] += fastevents->vweight[i]*Scale;
        
        Sys[erecbin] += getSystUncertainty(i,12)*Scale*SysScale;

        // fill signal
        if ( TMath::Abs(fastevents->vmode[i])<30){
          hErec[3]->Fill(fastevents->vfqenue[i],fastevents->vweight[i]);
          if (TMath::Abs(fastevents->vnutype[i])!=12){
            hErec[8]->Fill(fastevents->vfqenue[i],fastevents->vweight[i]);
          }
          else if (fastevents->vmode[i]==1){
            hErec[6]->Fill(fastevents->vfqenue[i],fastevents->vweight[i]);
          }
          else{
            hErec[7]->Fill(fastevents->vfqenue[i],fastevents->vweight[i]);           
          }
        }
        else{
          hErec[4]->Fill(fastevents->vfqenue[i],fastevents->vweight[i]);
          hErec[9]->Fill(fastevents->vfqenue[i],fastevents->vweight[i]);
        }
      }
    }
  }

  // add up figure of merit in each bin
  return calcFOM(Pow,Nev,Sys,nn);
 */


  // the new way:
  

  // summarize?
  if (flgPrintSummary){
    if (iplt==1){
      plots1 = new summaryPlots("nue");
      plots1->setLargeArray(fastevents);
      plots1->Init();
    }
    if (iplt==2){
      plots2 = new summaryPlots("nue");
      plots2->setLargeArray(fastevents);
      plots2->Init();
    }
  }

  // towall must be smaller than wall
  if (towallmin<wallmin){
    return 0.;
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

  // loop over events
  for (int i=0; i<nevents; i++){

     // see if event passes numu cuts
     int ipass = fastevents->vpassnue[i];

     // if it passes, add to spectrum
     if (ipass){

      // ..as long as it passes the FV cuts
      if ((fastevents->vfqwall[i] >= wallmin)&&(fastevents->vfqtowall[i]>=towallmin)){
       
        // get erec bin of this event
        int erecbin = hErec[0]->FindBin(fastevents->vfqenue[i]);

        // count the overflow bin
        if ((erecbin<0) && (fastevents->vfqenue[i]) > hErec[0]->GetBinLowEdge(hErec[0]->GetNbinsX())){
          erecbin = hErec[0]->GetNbinsX() + 1;
        }
       
        // get oscillation power (derivative w.r.t. oscillation parameter)
        float opow = getOscPowerFast(12,i,oscpar)*Scale;
        Pow[erecbin] += opow;
       
        // keep track of all events
        float ww = getEventWeight(i)*Scale;
        Nev[erecbin] += ww;
       
        // keep track of systematic uncertainty
        float syserr =  getSystUncertainty(i,12)*Scale*SysScale;
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

//////////////////////////////////////////
// get FOM using some erec binning 
float optimusPrime::calcFOMSpectrumNuMu(float towallmin, float wallmin, int oscpar, int iplt){

  // summarize?
  if (flgPrintSummary){
    if (iplt==1){
      plots1 = new summaryPlots("numuplt1");
      plots1->setLargeArray(fastevents);
      plots1->Init();
    }
    if (iplt==2){
      plots2 = new summaryPlots("numuplt2");
      plots2->setLargeArray(fastevents);
      plots2->Init();
    }
  }

  // towall must be smaller than wall
  if (towallmin<wallmin){
    return 0.;
  }

  // reset values
  Nevents = 0.;
  NB = 0.; 
  NS = 0.; 
  Power = 0.;
  Syst = 0.;

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

  // loop over events
  for (int i=0; i<nevents; i++){

     // see if event passes numu cuts
//     int ipass = passNuMuCuts(i);
     int ipass = fastevents->vpassnumu[i];

     // if it passes, add to spectrum
     if (ipass){

      // ..as long as it passes the FV cuts
      if ((fastevents->vfqwall[i] >= wallmin)&&(fastevents->vfqtowall[i]>=towallmin)){
       
        // get erec bin of this event
        int erecbin = hErec[0]->FindBin(fastevents->vfqenumu[i]);

        // count the overflow bin
        if ((erecbin<0) && (fastevents->vfqenumu[i]) > hErec[0]->GetBinLowEdge(hErec[0]->GetNbinsX())){
          erecbin = hErec[0]->GetNbinsX() + 1;
        }
        
        float opow = getOscPowerFast(14,i,oscpar)*Scale;
        Pow[erecbin] += opow;
        
        Nev[erecbin] += getEventWeight(i)*Scale;
        
        Sys[erecbin] += getSystUncertainty(i,14)*Scale*SysScale;

        if (flgPrintSummary){
          if (iplt==1) plots1->fillAllFromArray(i,0.,0.);
          if (iplt==2) plots2->fillAllFromArray(i,0.,0.);
        }

      }
    }
  }

  // add up figure of merit in each bin
  return calcFOM(Pow,Nev,Sys,nn);

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
// compare tow sets of FV cuts
void optimusPrime::compareFOM(float tw1, float w1, float tw2, float w2, int oscpar, int flgnumu){
  
  // make sure to record histograms
  flgPrintSummary = 1;

  // canvas setup
  multiPad = new TCanvas("multiPad","multiPad",700,800);
  multiPad->Divide(2,2);

  // colerz
  int colerz[5] = {4,7,2,15,6};

  // do the calculations and fill plots1
  float fom1 = 0.;

  if (flgnumu) fom1 = calcFOMSpectrumNuMu(tw1,w1,oscpar, 1);
  else{
    fom1 = calcFOMSpectrumNuE(tw1,w1,oscpar, 1);
  }

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
  if (flgnumu) fom2 = calcFOMSpectrumNuMu(tw2,w2,oscpar, 2);
  else{
    fom2 = calcFOMSpectrumNuE(tw2,w2,oscpar,2);
  }

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
void optimusPrime::compareCuts(float tw1, float w1, float tw2, float w2, int oscpar, int flgnumu){
  
  // make sure to record histograms
  flgPrintSummary = 1;

  // canvas setup
  multiPad = new TCanvas("multiPad","multiPad",700,800);
  multiPad->Divide(2,3);

  // colerz
  int colerz[5] = {4,7,2,15,6};

  // do the calculations and fill plots1
  float fom1 = 0.;

  if (flgnumu) fom1 = calcFOMSpectrumNuMu(tw1,w1,oscpar, 1);
  else{
    fom1 = calcFOMSpectrumNuE(tw1,w1,oscpar);
  }

  // do the calculations and fill plots2
  float fom2 = 0.;
  if (flgnumu) fom2 = calcFOMSpectrumNuMu(tw2,w2,oscpar, 2);
  else{
    fom2 = calcFOMSpectrumNuE(tw2,w2,oscpar);
  }

  // draw that shizz
  // set colors (loop of # of catagories in summaryPlots::GetCatagory() )

  if (flgnumu){ 
    multiPad->cd(1);
    plots1->pltEnuMu->Draw();
    plots1->pltEnuMu->SetLineWidth(3);
    plots2->pltEnuMu->SetLineWidth(3);
    plots2->pltEnuMu->SetLineStyle(2);
    plots1->pltEnuMu->Draw("h");
    plots2->pltEnuMu->Draw("same");
    double ymin = 0.;
    double ymax = plots1->pltEnuMu->GetMaximum();
    for (int icat=0; icat<5; icat++){
      plots1->pltEnuMuCat[icat]->SetLineColor(colerz[icat]);
      plots2->pltEnuMuCat[icat]->SetLineColor(colerz[icat]);
      plots2->pltEnuMuCat[icat]->SetLineWidth(3);
      plots1->pltEnuMuCat[icat]->SetLineWidth(3);
      plots2->pltEnuMuCat[icat]->SetLineStyle(2);
      multiPad->cd(icat+2);
//      plots1->pltEnuMuCat[icat]->SetMinimum(ymin);
//      plots1->pltEnuMuCat[icat]->SetMaximum(ymax);
      plots1->pltEnuMuCat[icat]->Draw("h");
      plots2->pltEnuMuCat[icat]->Draw("sameh");
    }
  }
  else{
    multiPad->cd(1);
    plots1->pltEnuE->Draw();
    plots1->pltEnuE->SetLineWidth(3);
    plots2->pltEnuE->SetLineWidth(3);
    plots2->pltEnuE->SetLineStyle(2);
    plots1->pltEnuE->Draw("h");
    plots2->pltEnuE->Draw("same");
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
//      plots1->pltEnuECat[icat]->SetMinimum(ymin);
//      plots1->pltEnuECat[icat]->SetMaximum(ymax);
      plots2->pltEnuECat[icat]->Draw("sameh");
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
// calculate FOM from arrays instead of histograms
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
      fom_thisbin = (pow[i]*pow[i])/((sys[i]*sys[i])+nev[i]);
    }
   
    if (flgPrintSummary){
      hErec[0]->SetBinContent(i,(pow[i]*pow[i]));
      hErec[1]->SetBinContent(i,nev[i]);
      hErec[2]->SetBinContent(i,sys[i]);
      hErec[3]->SetBinContent(i,fom_thisbin);
      hErec[4]->SetBinContent(i,pow[i]);
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
 }

 //
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

//   float lmom = fastevents->vfqenumu[i];
//   if (nutype!=14) lmom = fastevents->vfqenue[i];
   float sys = uncertaintyCalculator->getTotalUncertainty(fastevents->vwallv[i],
                                                          fastevents->vfqwall[i],
                                                          fastevents->vfqtowall[i],
                                                          fastevents->vfqenumu[i],
                                                          fastevents->vmode[i],
                                                          wronglepton);
//   if (sys>0.5)  cout<<"event "<<i<<" has sys: "<<sys<<endl;

   return sys*getEventWeight(i);                                                       
}



/////////////////////////////////////////////////////////////////
// use Toy MC to calculate systematics, instead of relying on maps
/*
float optimusPrime::calcFOMToyMC(float towallcut,
                                 float wallcut,
                                 int   oscpar,
                                 int flgselection,
                                 int nmcmcpts){

  // summarize?
  #ifdef PRINTSUMMARY
  plots = new summaryPlots("toyplots");
  plots->setLargeArray(fastevents);
  plots->InitToys(hErec[0]);
  #endif

  // reset histos
  hErec[0]->Reset();
  hErec[1]->Reset();
  hErec[2]->Reset();

  // set up shape parameters and classes to apply them
  atmFitPars* atmpars = new atmFitPars(cardFileName.Data()); //< create parameter containter 
  TChain* chmcmcpars = new TChain("MCMCpath");
  chmcmcpars->Add(mcmcParFileName.Data());
  mcmcReader* mcmcpars = new mcmcReader(chmcmcpars);
  modifier = new mcmcApply(atmpars, mcmcpars); //< this will perform the necessary modifications
  atmpars->resetDefaults();

  // make arrays to store info
  const int nbins = hErec[0]->GetNbinsX()+1;
  const int npoints = nmcmcpts;
  float A[npoints]; 
  float Powsq[nbins][npoints];
  float Nev[nbins][npoints];
  float Sys[nbins][npoints];

  // set initial array values
  for (int ibin=0; ibin<nbins; ibin++){
    for (int ipt=0; ipt<nmcmcpts; ipt++){
      Powsq[ibin][ipt]=0.;
      Nev[ibin][ipt]=0.;
      Sys[ibin][ipt]=0.;
    }
  }

  // get the selection of mcmc points
  vector<int> mcmcpoints;
  TRandom2* rando = new TRandom2(nmcmcpts);
  for (int ipt=0; ipt<nmcmcpts; ipt++){
    mcmcpoints.push_back( rando->Integer(chmcmcpars->GetEntries()) ); 
  }
  // sort that list
  std::sort(mcmcpoints.begin(),mcmcpoints.end());

  // do the toy MC
 
  // loop over the mcmc poins
  for (int ipt=0; ipt<nmcmcpts; ipt++){

    // get shape parameters
    cout<<"getting mcmc point: "<<mcmcpoints.at(ipt)<<endl;
    chmcmcpars->GetEntry(mcmcpoints.at(ipt));

    // set par containter values
    if (ipt>0) modifier->setFromMCMC();

    // loop over MC events 
    for (int iev=0; iev<nevents; iev++){

      // passes FV cuts?
      if ((fastevents->vfqtowall[iev]<towallcut)||(fastevents->vfqwall[iev]<wallcut)) continue;

      // apply the cuts to the modified event
//      int ipass = applyCutsToModifiedEvent(iev);
      int ipass = modifier->applyCutsToModifiedEvent(iev,fastevents);

      // is this the event you are looking for?
      if (ipass!=flgselection) continue;
      
      // get more info
      float enu=0.;
      float oscpow=0.;
      // for electron 
      if (flgselection==1){
        enu = fastevents->vfqenue[iev];
        oscpow = getOscPowerFast(12,iev,oscpar);
      }
      // for muon
      else if (flgselection==2){
        enu = fastevents->vfqenumu[iev];
        oscpow = getOscPowerFast(14,iev,oscpar);
      }

      // bin in energy
      int enubin = hErec[0]->FindBin(enu);
      hErec[0]->Fill(enu,fastevents->vweight[iev]);

      #ifdef PRINTSUMMARY
      plots->fillAllFromArray(iev,oscpow,0.);
      plots->pltToySpectrum[ipt]->Fill(enu,fastevents->vweight[iev]);
      plots->pltToyPower[ipt]->Fill(enu,oscpow*oscpow);
      #endif

      // fill arrays
      Powsq[enubin][ipt] += oscpow*oscpow;
      Nev[enubin][ipt] += fastevents->vweight[iev];

    } //< end loop over events

    // now calculate A (the curvature of likelihood for this selection and assumption of systematics)
//    for (int ibin=0; ibin<nbins; ibin++){
//      if (Nev[ibin][ipt]>0.) A[ipt] += Pow[ibin][ipt]/Nev[ibin][ipt];
//    }

  } //< end loop over mcmc points in toy MC


//  float meanA = arraymean(A,npoints);
//  cout<<"mean A: "<<meanA<<endl;
//  float varA  = arrayvar(A,npoints,meanA);
//  cout<<"var A: "<<varA<<endl;
//  hCurve = new TH1D("hA","hA",20,meanA-4*TMath::Sqrt(varA),
//                                 meanA+4*TMath::Sqrt(varA));
//  for (int ipt=0; ipt<npoints; ipt++){
//    hCurve->Fill(A[ipt]);
//  }

   // loop over bins to calculate FOM
   float FOM = 0.;
   for (int ibin=0; ibin<=nbins; ibin++){
      
      // get bin mean
      float nominal_bin_content = Nev[ibin][0];
      cout<<"nominal contnet: "<<nominal_bin_content<<endl;
      float mean_bin_content = arraymean(Nev[ibin],npoints);
      cout<<"mean contnet: "<<mean_bin_content<<endl;
      float nominal_bin_power = Powsq[ibin][0];
      cout<<"nominal bin power: "<<nominal_bin_power<<endl;
      float nominal_bin_syst  = Sys[ibin][0];
      cout<<"nominal bin syst: "<<nominal_bin_syst<<endl;
      float bin_content_variance = arrayvar(Nev[ibin],npoints,mean_bin_content);
      cout<<"bin content var : "<<bin_content_variance<<endl;
      float bin_content_shift = mean_bin_content - nominal_bin_content;
      float unc_total = bin_content_variance + nominal_bin_syst*nominal_bin_syst + ((bin_content_shift)*(bin_content_shift));
      if (nominal_bin_content>0.) FOM += (nominal_bin_power)/(nominal_bin_content + unc_total);

//      hErec[0]->SetBinContent(ibin,mean_bin_content);
      hErec[0]->SetBinContent(ibin,nominal_bin_content);
      hErec[0]->SetBinError(ibin,TMath::Sqrt(bin_content_variance));
//      hErec[1]->SetBinContent(ibin,nominal_bin_content);
      hErec[1]->SetBinContent(ibin,nominal_bin_power);


   }
 

  // calculate mean and covariance
//  float meanPow[nbins];
//  float meanNev[nbins];
//  float cov[nbins][nbins];

  // init to zero
//  for (int ibin=0; ibin<nn; ibin++){
//    mean[ibin] = 0.;
//    for (int jbin=0; jbin<nn; jbin++){
//      cov[ibin][jbin] = 0.;
//    }
//  }

//  for (int ibin=0; ibin<nn; ibin++){
//    for (int ipt=0; ipt<nmcmcpts; ipt++){
//       mean[ibin] +=  

  return FOM;
}
*/


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

  /*
  cutPars.fqmommu = fastevents->vfqmumom[i]; 
  cutPars.fqmome = fastevents->vfqemom[i];
  cutPars.fqpid = fastevents->vattribute[i][indexPIDPar];
  cutPars.fqpi0par =fastevents->vattribute[i][indexPi0Par];
  cutPars.fqpippar = fastevents->vattribute[i][indexPiPPar];
  cutPars.fqenumu = fastevents->vfqenumu[i];
  cutPars.fqenue = fastevents->vfqenue[i];
  cutPars.fqrcpar = fastevents->vfqrcpar[i];
  cutPars.nhitac = fastevents->vnhitac[i];
  cutPars.fqnsubev = fastevents->vfqnsubev[i];
  cutPars.fqnring = fastevents->vfqnring[i];

 
  int ipass = selectNuMu(cutPars);
  */

//  int ipass = selectNuMu( fastevents->vnhitac[i],
//                          fastevents->vfqnsubev[i],
//                          fastevents->vfqenumu[i],
//                          fastevents->vfqemom[i],
//                          fastevents->vfqmumom[i],
//                          fastevents->vfqpid[i],
//                          fastevents->vfqnring[i] )


  int ipass = fastevents->vpassnumu[i];
  return ipass;

}

/////////////////////////////////////////
// what about nu-e cuts?
int optimusPrime::passNuECuts(int i){

  /*
  cutPars.fqmommu = fastevents->vfqmumom[i]; 
  cutPars.fqmome = fastevents->vfqemom[i];
  cutPars.fqpid = fastevents->vattribute[i][indexPIDPar];
  cutPars.fqpi0par =fastevents->vattribute[i][indexPi0Par];
  cutPars.fqpippar = fastevents->vattribute[i][indexPiPPar];
  cutPars.fqenumu = fastevents->vfqenumu[i];
  cutPars.fqenue = fastevents->vfqenue[i];
  cutPars.fqrcpar = fastevents->vfqrcpar[i];
  cutPars.fqnring = fastevents->vfqnring[i];
  cutPars.nhitac = fastevents->vnhitac[i];
  cutPars.fqnsubev = fastevents->vfqnsubev[i];


  int ipass = selectNuE(cutPars);

//  int ipass = selectNuE(  fastevents->vnhitac[i],
//                          fastevents->vfqnsubev[i],
//                          fastevents->vfqenumu[i],
//                          fastevents->vfqemom[i],
//                          fastevents->vfqpid[i],
//                          fastevents->vfqnring[i],
//                          fastevents->vfqpi0par[i]);
*/                          
  int ipass = fastevents->vpassnue[i];
  return ipass;

}


/////////////////////////////////////////
// get total FOM using single bin
float optimusPrime::calcNuMuFOM(float towallmin, float wallmin, int oscpar){

  if (towallmin<wallmin) return 0.;

  Nevents = 0.;
  NB = 0.; 
  NS = 0.; 
  Power = 0.;
  Syst = 0.;

  for (int i=0; i<nevents; i++){
     int ipass = passNuMuCuts(i);
//     int ipass = selectNuMu( fastevents->vnhitac[i],
//                             fastevents->vfqnsubev[i],
//                             fastevents->vfqenumu[i],
//                             fastevents->vfqemom[i],
//                             fastevents->vfqmumom[i],
//                             fastevents->vfqpid[i],
//                             fastevents->vfqnring[i] );
     if (ipass){
      if ((fastevents->vfqwall[i] > wallmin)&&(fastevents->vfqtowall[i]>towallmin)){
        if (fastevents->vfqenumu[i]>10000.) continue;
        Power += getOscPowerFast(14,i,oscpar);
        Nevents += fastevents->vweight[i];
        Syst += getSystUncertainty(i);
        if ((fastevents->vmode[i]!=1 || fastevents->vnutype[i]!=14)){
          NB+=fastevents->vweight[i];
        }
        else{
          NS+=fastevents->vweight[i];
        }
      }
    }
  }

  Power*=Scale;
  Nevents*=Scale;
  NS*=Scale;
  NB*=Scale;
  Syst*=Scale;

  float fom = (Power*Power)/((Nevents+(Syst*Syst)));
  if (FOMType==1) fom = Nevents;
  if (FOMType==2) fom = Syst*Syst;
  if (FOMType==3) fom = Power*Power/(Nevents*Nevents);
  if (FOMType==4) fom = Power*Power;
  if (FOMType==5) fom = NS;
  if (FOMType==6) fom = NB;
  if (FOMType==7) fom = NS/NB;
  if (FOMType==8) fom = NS/TMath::Sqrt(NB+NS);
  if (FOMType==9) fom = (NS*NS)/Nevents;
  if (FOMType==10) fom = (Syst)/Nevents;



//  float scale = (float)nevents/(float)chmc->GetEntries();
//  fom*=1./scale;
  
  cout<<"FOM: "<<fom<<endl;
  cout<<"N: "<<Nevents<<endl;
  cout<<"S: "<<Syst<<endl;
  cout<<"P: "<<Power<<endl;

  return fom;

}

#endif







