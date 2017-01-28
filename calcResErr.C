///////////////////////////////////////////////////////////////////////////
// Functions to calculate "residual" uncertainty between two mc histograms
#ifndef CALCRES_C
#define CALCRES_C

#include "TH1D.h"

#include <iostream>
#include "TMath.h"

using namespace std;

void cumulatebackward(TH1D* hh){
  
  const int N = hh->GetNbinsX();
  double cc[N];

  int lastbin = hh->GetNbinsX();
  for (int ibin=lastbin; ibin<=1; ibin--){
    cc[ibin] = hh->Integral(1,lastbin); 
  }

  return;
}

double integrateError(TH1D* hh, int binlow, int binhigh){
  
  double sumerr = 0.;
  for (int i=binlow; i<=binhigh; i++){
     double binerr = hh->GetBinError(i);
     sumerr+=(binerr*binerr);
  }
  
  return TMath::Sqrt(sumerr);
}

void cumulateforward(TH1D* hh){
  
  const int N = hh->GetNbinsX()+1;
  double cc[N];
  double ccerr[N];

  int lastbin = hh->GetNbinsX();
  for (int ibin=0; ibin<=lastbin; ibin++){
    cc[ibin] = hh->Integral(ibin,lastbin); 
    ccerr[ibin] = integrateError(hh,ibin,lastbin);
  }
  
  for (int ibin=0; ibin<=lastbin; ibin++){
    hh->SetBinContent(ibin,cc[ibin]);
    hh->SetBinError(ibin,ccerr[ibin]);
  }
  

  return;
}

int getMostSigBin(TH1D* hh, double baseline=0.0){

  int nn = hh->GetNbinsX();

  int isigbin = -1;
  double sigval = -1.;
  double sigbest = -1;
  for (int ibin=1; ibin<nn; ibin++){
    if (hh->GetBinError(ibin)<=0.) continue;
    sigval = (hh->GetBinContent(ibin) - baseline)/(hh->GetBinError(ibin));
    if (TMath::Abs(sigval)>sigbest){
      sigbest = TMath::Abs(sigval);
      isigbin = ibin;
    }
  }


  return isigbin;
}


void calcCumDiff(TH1D* hmc, TH1D* hdata){

   cumulateforward(hmc);
   cumulateforward(hdata);
   hmc->Add(hdata,-1.);
   for (int ibin=1; ibin<=hmc->GetNbinsX(); ibin++){
     double binc = hmc->GetBinContent(ibin);
     if (binc<0.){
       hmc->SetBinContent(ibin,TMath::Abs(binc)); 
     }
   }

   return;
}

void calcCumFracDiff(TH1D* hmc, TH1D* hdata){

   cumulateforward(hmc);
   cumulateforward(hdata);
   TH1D* hclone = (TH1D*)hmc->Clone("hclonemc");
   hmc->Add(hdata,-1.);
   hmc->Divide(hclone); 

   for (int ibin=1; ibin<=hmc->GetNbinsX(); ibin++){
     double binc = hmc->GetBinContent(ibin);
     if (binc<0.){
       hmc->SetBinContent(ibin,TMath::Abs(binc)); 
     }
   }

   return;
}


double calcFracDiff(TH1D* hmc, TH1D* hdata){
  
  TH1D* hclonem = (TH1D*)hmc->Clone("hclonem");
  TH1D* hcloned = (TH1D*)hdata->Clone("hcloned");

  hclonem->Add(hcloned,-1.);
  hclonem->Divide(hmc);

  double wsum = 0.;
  double norm = 0.;
  for (int ibin=1; ibin<=hmc->GetNbinsX(); ibin++){
    wsum+=(hmc->GetBinContent(ibin))*TMath::Abs(hclonem->GetBinContent(ibin));
    norm+= hmc->GetBinContent(ibin);
    hmc->SetBinContent(ibin,hmc->GetBinContent(ibin)*TMath::Abs(hclonem->GetBinContent(ibin)));
  }

  wsum/=norm;

//  hmc=(TH1D*)hclonem->Clone("hfinal");
  return wsum;
 }

// integrate method
double calcCumFracDiffI(TH1D* hmc, TH1D* hdata){
   TH1D* hclone0 = (TH1D*)hmc->Clone("hclonemc0");
   TH1D* hclonedata0 = (TH1D*)hdata->Clone("hclonedata0");

   cumulateforward(hclone0);
   cumulateforward(hclonedata0);
   TH1D* hclone1 = (TH1D*)hclone0->Clone("hclone1");
   hclone0->Add(hclonedata0,-1.);
   hclone0->Divide(hclone1); 

   double binsum = 0.;
   double norm = 1./hmc->Integral();
   for (int ibin=1; ibin<=hmc->GetNbinsX(); ibin++){
     double binc = hclone0->GetBinContent(ibin);
     binc *= hmc->GetBinContent(ibin)*norm;
     if (binc<0.){
       binc = TMath::Abs(binc); 
     }
     binsum += binc;
 //    hclone0->SetBinContent(ibin,binc);
  //   hmc->SetBinError(ibin,0.);
   }

   hclone0->Delete();
   hclone1->Delete();
   hclonedata0->Delete();
   return binsum;
}


double calcResErrorI(TH1D* hmc, TH1D* hdata){
  
  // clone 
  TH1D* hclone = (TH1D*)hmc->Clone("hclone");
  TH1D* hcloned = (TH1D*)hdata->Clone("hcloned");

  // get abosolute cumulative difference
  calcCumFracDiffI(hclone,hcloned);
 
  // return integral 
  double avg = hclone->Integral();

  hclone->Delete();
  hcloned->Delete();

  return avg; 

}



/*
double calcResError(TH1D* hmc, TH1D* hdata){
  
  // clone 
  TH1D* hclone = (TH1D*)hmc->Clone("hclone");
  TH1D* hcloned = (TH1D*)hdata->Clone("hcloned");

  // get abosolute cumulative difference
  calcCumFracDiff(hclone,hcloned);
 
  // return max value
  double max = hclone->GetMaximum(); 

  hclone->Delete();
  hcloned->Delete();

  return max; 

}
*/


double calcResErrorFrac(TH1D* hmc, TH1D* hdata){
  
  // clone 
  TH1D* hclone = (TH1D*)hmc->Clone("hclone");
  TH1D* hcloned = (TH1D*)hdata->Clone("hcloned");

  // get abosolute cumulative difference
  calcCumDiff(hclone,hcloned);
 
  // return max value
  double max = hclone->GetMaximum(); 

  hclone->Delete();
  hcloned->Delete();

  return max; 

}

void absVal(TH1D* hh){

  int nn = hh->GetNbinsX();
  for (int i=1; i<=nn; i++){
    double binc = hh->GetBinContent(i);
    if (hh->GetBinContent(i)<0.){
      hh->SetBinContent(i,TMath::Abs(binc));
    }
  }

  return;
  
}

double calcResErrorFrac2(TH1D* hmc, TH1D* hdata){
  
  // clone 
  TH1D* hclone = (TH1D*)hmc->Clone("hclone");
  TH1D* hcloned = (TH1D*)hdata->Clone("hcloned");

//  hclone->Add(hcloned,-1.);
//  hclone->Divide(hmc);
//  absVal(hclone);
 

  hmc->Add(hdata,-1.);
  hmc->Divide(hclone);
  absVal(hmc); 
  hclone->Delete();
  hcloned->Delete();

  return 1.0; 

}


double calcResErrorFrac3(TH1D* hmc, TH1D* hdata){
  
  // clone 
  TH1D* hclone = (TH1D*)hmc->Clone("hclone");
  TH1D* hcloned = (TH1D*)hdata->Clone("hcloned");

  hmc->Add(hdata,-1.);
//  hmc->Scale(1./hclone->Integral());
//  hmc->Divide(hclone);
 //j: absVal(hmc); 
  hclone->Delete();
  hcloned->Delete(); 


  return 1.0; 

}



void getSignif(TH1D* hmc, TH1D* hdata){
  hmc->Sumw2(); 
  hdata->Sumw2(); 

  TH1D* hclone = (TH1D*)hmc->Clone("hclone");
  TH1D* hcloned = (TH1D*)hdata->Clone("hcloned");

  hclone->Divide(hdata);
  
  for (int ibin=1; ibin<=hmc->GetNbinsX(); ibin++){
    double diff = hclone->GetBinContent(ibin) -1.;
    double unc = hclone->GetBinError(ibin);
    cout<<diff<<endl;
    cout<<unc<<endl;
    cout<<"--------------"<<endl;
    if (unc>0.) hmc->SetBinContent(ibin,TMath::Abs(diff)/unc);
  }

  return ;
}


double calcResError(TH1D* hmc, TH1D* hdata){
  return calcFracDiff(hmc,hdata);
//  return calcCumFracDiffI(hmc,hdata);
}




void applyResError(TH1D* hmc, TH1D* hdata){

  double fracerr = calcResError(hmc,hdata)/hmc->Integral();

  for (int ibin=1; ibin<hmc->GetNbinsX(); ibin++){
     double binerr = hmc->GetBinError(ibin); 
     double binc = hmc->GetBinContent(ibin); 
     double reserr = binc*fracerr;
     double toterr = TMath::Sqrt(binerr*binerr + reserr*reserr);
     hmc->SetBinError(ibin,toterr);
  }

  return;
}

#endif
