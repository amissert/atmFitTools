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
// if (TMath::Abs(fastevents->vmode[ientry])>=30) return 0.;
 if (TMath::Abs(fastevents->vmode[ientry])!=1) return 0.;
 
 // other nu do not contribute
 if (TMath::Abs(fastevents->vnutype[ientry] != nutype)) return 0.;

 // outside FV do not contribute
 if (fastevents->vwallv[ientry] <=0.) return 0.;

 // return power multiplied by weight
 float oscpow = fastevents->vweight[ientry]*fastevents->voscpower[ientry][oscpar]; 

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

void optimusPrime::calcFOMMap(float towallmax, float wallmax,int oscpar,int npts){

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
   else{ value = calcNuMuFOMSpectrum(towallcut,wallcut,oscpar);}
   if (value>fommax){
     fommax = value;
     maxbin = ibin;
   }
   hFV->SetBinContent(ibin,value);
 }
 
// float frange = 0.1*(fmax-fmin);
// hFV->GetZaxis()->SetRangeUser(fmin-frange,fmax+frange);
// hFV->GetZaxis()->SetRange(fmin-frange,fmax+frange);
// hFV->SetAxisRange(fmin-frange,fmax+frange,"Z");
// hFV->SetMaximum(fmax+frange);
// hFV->SetMinimum(fmax+frange);
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


//////////////////////////////////////////
// get FOM using a spectrum
float optimusPrime::calcNuMuFOMSpectrum(float towallmin, float wallmin, int oscpar){

  if (towallmin<wallmin){
    return 0.;
  }

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
     int ipass = passNuMuCuts(i);

     // if it passes, add to spectrum
     if (ipass){

      // ..as long as it passes the FV cuts
      if ((fastevents->vfqwall[i] >= wallmin)&&(fastevents->vfqtowall[i]>=towallmin)){
       
        // get erec bin of this event
//        if (fastevents->vfqenumu[i]>10000) continue;
        int erecbin = hErec[0]->FindBin(fastevents->vfqenumu[i]);
         
        // add to power
//        hErec[0]->Fill(fastevents->vfqenumu[i], getOscPowerFast(14,i,oscpar));
        Pow[erecbin] += getOscPowerFast(14,i,oscpar)*Scale;
        
        // add to nevents
//        hErec[1]->Fill(fastevents->vfqenumu[i],fastevents->vweight[i]);
        Nev[erecbin] += fastevents->vweight[i]*Scale;
        
        // add to systematics
//        if (hErec[2]->Fill(fastevents->vfqenumu[i], getSystUncertainty(i))<0){
//          cout<<"bad EREC: "<<fastevents->vfqenumu[i]<<endl;
//        }
        Sys[erecbin] += getSystUncertainty(i)*Scale*SysScale;
//        Syst+=getSystUncertainty(i);

        // fill signal
        if ( TMath::Abs(fastevents->vmode[i])==1){
          hErec[3]->Fill(fastevents->vfqenumu[i],fastevents->vweight[i]);
        }
        // fill background
        else{
          hErec[4]->Fill(fastevents->vfqenumu[i],fastevents->vweight[i]); 
        }
      }
    }
  }

  // add up figure of merit in each bin
//  return calcFOMSpectrum(Pow,Nev,Sys,nn);
  return TMath::Sqrt(calcFOMSpectrum(Pow,Nev,Sys,nn));
//  cout<<"sumsyst: "<<Syst<<endl;
//  return calcFOMSpectrum();
}

/////////////////////////////////////////////////////
// compare tow sets of FV cuts
void optimusPrime::compareNuMuCuts(float tw1, float w1, float tw2, float w2, int oscpar){
  
  // canvas setup
  multiPad = new TCanvas("multiPad","multiPad",700,800);
  multiPad->Divide(2,3);

  float fom1 = calcNuMuFOMSpectrum(tw1,w1,oscpar);
  for (int ih=0; ih<6; ih++){
    multiPad->cd(ih+1);
    TString hname = Form("hs%d",ih);
    hSummary[ih] = (TH1D*)hErec[ih]->Clone(hname.Data());  
    hSummary[ih]->SetLineColor(kRed);  
    hSummary[ih]->SetLineWidth(3);  
    hSummary[ih]->Draw("h");
  }

  float fom2 = calcNuMuFOMSpectrum(tw2,w2,oscpar);
    for (int ih=0; ih<6; ih++){
    multiPad->cd(ih+1);
    hErec[ih]->SetLineWidth(3);  
    hErec[ih]->Draw("sameh");
  } 


  //
  return;
}

///////////////////////////////////////////////////
// calculate FOM from arrays instead of histograms
float optimusPrime::calcFOMSpectrum(float* pow, float* nev, float* sys, int nbin){

 float fom = 0.; 
 float S  = 0.;
 float N = 0.;
 float P = 0.;

 for (int i=0; i<=nbin; i++){
    float fombin = 0.;
    S+=sys[i];
    N+=nev[i];
    P+=pow[i];
    if (FOMType==0) fombin = (pow[i]*pow[i])/((sys[i]*sys[i])+nev[i]);
    if (FOMType==1) fombin = (pow[i]);
    if (FOMType==2) fombin = (sys[i]);  
    if (FOMType==3) fombin = (nev[i]);  
    if ((sys[i]+nev[i])>0.){
      hErec[5]->SetBinContent(i,fombin);
      fom+=fombin;
    }
    else{
      hErec[5]->SetBinContent(i,0.);
    }
    hErec[0]->SetBinContent(i,TMath::Abs(pow[i]));
    hErec[1]->SetBinContent(i,nev[i]);
    hErec[2]->SetBinContent(i,sys[i]);
 }
 
 cout<<"P: "<<P<<endl;
 cout<<"N: "<<N<<endl;
 cout<<"S: "<<S<<endl;
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
     continue;
   }
   // fill 
   float ww = fastevents->vweight[iev];
   float sys = getSystUncertainty(iev); 
   float pow = getOscPowerFast(nutype,iev,oscpar);
   float enurc = fastevents->vfqenumu[iev];
   float enuv = fastevents->vpmomv[iev];
   // use RC
//   hFVSummary[0]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],ww);
//   hFVSummary[1]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],pow);
//   hFVSummary[2]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],sys);
//   hFVSummary[4]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],enuv*ww);
//   hFVSummary[5]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],enurc*ww);
   // use true
   hFVSummary[0]->Fill(fastevents->vtowallv[iev],fastevents->vwallv[iev],ww);
   hFVSummary[1]->Fill(fastevents->vtowallv[iev],fastevents->vwallv[iev],pow);
   hFVSummary[2]->Fill(fastevents->vtowallv[iev],fastevents->vwallv[iev],sys);
   hFVSummary[4]->Fill(fastevents->vtowallv[iev],fastevents->vwallv[iev],enuv*ww);
   hFVSummary[5]->Fill(fastevents->vtowallv[iev],fastevents->vwallv[iev],enurc*ww);

   if (TMath::Abs(fastevents->vmode[iev])==1) hFVSummary[6]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],ww);
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
float optimusPrime::getSystUncertainty(int i){
   float sys = uncertaintyCalculator->getTotalUncertainty(fastevents->vwallv[i],
                                                          fastevents->vfqwall[i],
                                                          fastevents->vfqtowall[i],
                                                          fastevents->vfqenumu[i],
                                                          fastevents->vmode[i])*fastevents->vweight[i];
   return sys;                                                       
}


/////////////////////////////////////////
// use the filled hErec histos to calculate
// the FOM for a set of cuts
float optimusPrime::calcFOMSpectrum(){
  
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


///////////////////////////////////////////
// does event pass cuts?
int optimusPrime::passNuMuCuts(int i){

  int ipass = selectNuMu( fastevents->vnhitac[i],
                          fastevents->vfqnsubev[i],
                          fastevents->vfqenumu[i],
                          fastevents->vfqemom[i],
                          fastevents->vfqmumom[i],
                          fastevents->vfqpid[i],
                          fastevents->vfqnring[i] );
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







