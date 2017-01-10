#ifndef PRIME_CXX
#define PRIME_CXX


#include "optimusPrime.h"

//////////////////////////////////////////
// constructor
optimusPrime::optimusPrime(TChain* t2kmc, int nevts){

 chmc = t2kmc;
 mcevent = new fqProcessedEvent(chmc);
 nevents = nevts;
 FOMType = 0;

 if (nevts>chmc->GetEntries()){
   nevents = chmc->GetEntries();
   flgUseEventList = 1;
 }
 else{
  flgUseEventList = 1;
 }

 // we only need enu, wall, towall, and mode and wallv
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

 hFVAll = new TH2FV("hall",-1,30,0,800,30,0,800);
 hFVAvg = new TH2FV("havg",-1,30,0,800,30,0,800);

 uncertaintyCalculator = new moreUncertainties("/nfs/data41/t2k/amissert/atmos/head/atmFitTools/data/");

 Scale = 1.;
 SysScale = 1.;

 fillArray();

 // setup recon energy histo
 int nbins = 100;
 float emax = 5000;

 TH1D* hseed = uncertaintyCalculator->hERecUnc[0];
 hseed->SetStats(0);
 hseed->SetBit(TH1::kNoTitle,0);
 hErec[0] = (TH1D*)hseed->Clone("herec_power");
// hErec[0] = new TH1D("h","h",20,0,2000);
// hseed = hErec[0];
 hErec[0]->SetTitle("Oscillation Power");
 hErec[1] = (TH1D*)hseed->Clone("herec_N");
 hErec[1]->SetTitle("# of Events");
 hErec[2] = (TH1D*)hseed->Clone("herec_syst");
 hErec[2]->SetTitle("Uncertainty");
 hErec[3] = (TH1D*)hseed->Clone("herec_signal");
 hErec[3]->SetTitle("# Signal");
 hErec[4] = (TH1D*)hseed->Clone("herec_bg");
 hErec[4]->SetTitle("# BG");
 hErec[5] = (TH1D*)hseed->Clone("herec_fom");
 hErec[5]->SetTitle("F.O.M.");
 hErec[6] = (TH1D*)hseed->Clone("herec_ccqe");
 hErec[6]->SetTitle("CCQE");
 hErec[7] = (TH1D*)hseed->Clone("herec_ccnqe");
 hErec[7]->SetTitle("CCnQE");
 hErec[8] = (TH1D*)hseed->Clone("herec_ccwrong");
 hErec[8]->SetTitle("CCWrong");
 hErec[9] = (TH1D*)hseed->Clone("herec_nc");
 hErec[9]->SetTitle("NC");

 for (int ih=0; ih<6; ih++){
    hErec[ih]->Reset();
 }

 flgUseSpectrum = 1;
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

 // NC events do not contribute 
 if (TMath::Abs(fastevents->vmode[ientry])>=30) return 0.;
 
 // other nu do not contribute
 if (TMath::Abs(fastevents->vnutype[ientry]) != nutype) return 0.;

 // outside FV do not contribute
 if (fastevents->vwallv[ientry] <=0.) return 0.;

 // return power multiplied by weight
 float oscpow = fastevents->vweight[ientry]*fastevents->voscpower[ientry][oscpar]; 

// cout<<"oscpow: "<<oscpow<<endl;

 return oscpow;
}

void optimusPrime::fillFVHistoFast(){
//  hFVAll = new TH2FV("hall",-1,40,0,1000,40,0,1000);
//  hFVAvg = new TH2FV("havg",-1,40,0,1000,40,0,1000);
  hFVAll->Reset();
  hFVAvg->Reset();
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

 double fommax = -1.;
 int    maxbin = -1;

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
float optimusPrime::calcFOMSpectrumNuE(float towallmin, float wallmin, int oscpar){

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

}


/////////////////////////////////////////////
// draw stacked histogram of events
void optimusPrime::showBreakdown(){

  hs = new THStack("hs","");

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
// get FOM using a spectrum
float optimusPrime::calcFOMSpectrumNuMu(float towallmin, float wallmin, int oscpar){

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

  // nue or numu analysis?

  // loop over events
  for (int i=0; i<nevents; i++){

     // see if event passes numu cuts
     int ipass = passNuMuCuts(i);

     // if it passes, add to spectrum
     if (ipass){

      // ..as long as it passes the FV cuts
      if ((fastevents->vfqwall[i] >= wallmin)&&(fastevents->vfqtowall[i]>=towallmin)){
       
        // get erec bin of this event
        int erecbin = hErec[0]->FindBin(fastevents->vfqenumu[i]);
        
        Pow[erecbin] += getOscPowerFast(14,i,oscpar)*Scale;
        
        Nev[erecbin] += fastevents->vweight[i]*Scale;
        
        Sys[erecbin] += getSystUncertainty(i,14)*Scale*SysScale;

        // fill signal
        if ( TMath::Abs(fastevents->vmode[i])<30){
          hErec[3]->Fill(fastevents->vfqenumu[i],fastevents->vweight[i]);
        }
        else{
          hErec[4]->Fill(fastevents->vfqenumu[i],fastevents->vweight[i]);
        }
      }
    }
  }

  // add up figure of merit in each bin
  return calcFOM(Pow,Nev,Sys,nn);

}

/////////////////////////////////////////////////////
// compare tow sets of FV cuts
void optimusPrime::compareCuts(float tw1, float w1, float tw2, float w2, int oscpar,int flgnumu){;
  
  // canvas setup
  multiPad = new TCanvas("multiPad","multiPad",700,800);
  multiPad->Divide(2,3);

  float fom1 = 0.;
  if (flgnumu) fom1 = calcFOMSpectrumNuMu(tw1,w1,oscpar);
  else{
    fom1 = calcFOMSpectrumNuE(tw1,w1,oscpar);
  }

  for (int ih=0; ih<6; ih++){
    multiPad->cd(ih+1);
    TString hname = Form("hs%d",ih);
    hSummary[ih] = (TH1D*)hErec[ih]->Clone(hname.Data());  
    hSummary[ih]->SetLineColor(kRed);  
    hSummary[ih]->SetLineWidth(3);  
    hSummary[ih]->GetYaxis()->SetRangeUser(0,hSummary[ih]->GetMaximum());
    hSummary[ih]->Draw("h");
  }

  float fom2 = 0.;
  if (flgnumu) fom2 = calcFOMSpectrumNuMu(tw2,w2,oscpar);
  else{
    fom2 = calcFOMSpectrumNuE(tw2,w2,oscpar);
  }

  for (int ih=0; ih<6; ih++){
    multiPad->cd(ih+1);
    hErec[ih]->SetLineWidth(3);  
    hErec[ih]->Draw("sameh");
  } 


  //
  cout<<"FOM1: "<<fom1<<endl;
  cout<<"FOM2: "<<fom2<<endl;
  return;
}

///////////////////////////////////////////////////
// calculate FOM from arrays instead of histograms
float optimusPrime::calcFOM(float* pow, float* nev, float* sys, int nbin){

 float fom = 0.; 
 float S  = 0.;
 float N = 0.;
 float P = 0.;

 for (int i=0; i<nbin; i++){
    float fombin = 0.;
    S+=sys[i];
    N+=nev[i];
    P+=pow[i];
    if (FOMType==0){
      if ((sys[i]+nev[i])>0.){
        fombin = (pow[i]*pow[i])/((sys[i]*sys[i])+nev[i]);
//        fom+=fombin;
        hErec[5]->SetBinContent(i,fombin);
      }
      else{
        hErec[5]->SetBinContent(i,0.);
      }
    }
    if (FOMType==1) fombin = (pow[i]);
    if (FOMType==2) fombin = (sys[i]);  
    if (FOMType==3) fombin = (nev[i]);  
    fom+=fombin;
    hErec[0]->SetBinContent(i,TMath::Abs(pow[i]));
    hErec[1]->SetBinContent(i,nev[i]);
    hErec[2]->SetBinContent(i,sys[i]);
 }
 
 cout<<"total P: "<<P<<endl;
 cout<<"total N: "<<N<<endl;
 cout<<"total S: "<<S<<endl;
 cout<<"FOM: "<<fom<<endl;
 return fom;

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
   int wronglepton = 0;
   if (TMath::Abs(fastevents->vnutype[i])!=nutype) wronglepton = 1;
   float lmom = fastevents->vfqenumu[i];
   if (nutype!=14) lmom = fastevents->vfqenue[i];
   float sys = uncertaintyCalculator->getTotalUncertainty(fastevents->vwallv[i],
                                                          fastevents->vfqwall[i],
                                                          fastevents->vfqtowall[i],
                                                          fastevents->vfqenumu[i],
                                                          fastevents->vmode[i],
                                                          wronglepton);
//   if (sys>0.5)  cout<<"event "<<i<<" has sys: "<<sys<<endl;

   return sys*fastevents->vweight[i];                                                       
}


/////////////////////////////////////////
// use the filled hErec histos to calculate
// the FOM for a set of cuts
/*float optimusPrime::calcFOMSpectrum(){
  
  float fom = 0.;

  // loop over bins
  for (int ibin=1; ibin<=hErec[0]->GetNbinsX(); ibin++){
    float P = hErec[0]->GetBinContent(ibin);
    float N = hErec[1]->GetBinContent(ibin);
    float S = hErec[2]->GetBinContent(ibin);
    // add to F.O.M.
    if ((N+S)>0)fom += (P*P)/(N+(S*S));
  }

  //
  cout<<"P:   "<<hErec[0]->Integral()<<endl;
  cout<<"S:   "<<hErec[2]->Integral()<<endl;
  cout<<"N:   "<<hErec[1]->Integral()<<endl;
  cout<<"FOM: "<<fom<<endl;
  return fom;
}
*/


/////////////////////////////////////////////////////////////////
// use Toy MC to calculate systematics, instead of relying on maps
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

  cutPars.fqmommu = fastevents->vfqmumom[i]; 
  cutPars.fqmome = fastevents->vfqemom[i];
  cutPars.fqpid = fastevents->vattribute[i][indexPIDPar];
  cutPars.fqpi0par =fastevents->vattribute[i][indexPi0Par];
  cutPars.fqpippar = fastevents->vattribute[i][indexPiPPar];
  cutPars.fqenumu = fastevents->vfqenumu[i];
  cutPars.fqenue = fastevents->vfqenue[i];
  cutPars.fqrcpar = fastevents->vattribute[i][indexRCPar];
  cutPars.nhitac = fastevents->vnhitac[i];
  cutPars.fqnsubev = fastevents->vfqnsubev[i];

 
  int ipass = selectNuMu(cutPars);

//  int ipass = selectNuMu( fastevents->vnhitac[i],
//                          fastevents->vfqnsubev[i],
//                          fastevents->vfqenumu[i],
//                          fastevents->vfqemom[i],
//                          fastevents->vfqmumom[i],
//                          fastevents->vfqpid[i],
//                          fastevents->vfqnring[i] );
  return ipass;

}

/////////////////////////////////////////
// what about nu-e cuts?
int optimusPrime::passNuECuts(int i){

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
//                          
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







