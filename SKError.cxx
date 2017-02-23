#ifndef SKERROR_C
#define SKERROR_C

#include "SKError.h"



///////////////////////////////////////////////////////////////
float SKError::calcEff(int nclass, int ntoy, int totindex){

  // don't divide by zero
  if (NeventsTotal[nclass][totindex]==0) return 0.;
  if (NeventsTotal[nclass][ntoy]==0) return 0.;

  // get this efficiency  
  float eff = Nevents[nclass][ntoy]/NeventsTotal[nclass][ntoy];
  float eff0 = Nevents[nclass][totindex]/NeventsTotal[nclass][totindex];
  
  if (eff0 = 0.) return 0.;

  return ((eff/eff0) - 1.); 

}


///////////////////////////////////////////////////////////////
void SKError::drawScatter(int iclass, int jclass){

  const int nn = Ntoys;
  double X[nn];
  double Y[nn];

  for (int i=0; i<Ntoys; i++){
    X[i] = Nevents[iclass][i];
    Y[i] = Nevents[jclass][i];
  }

  gScat = new TGraph(nn,X,Y);

  gScat->Draw("ap");

  return;
}



///////////////////////////////////////////////////////////////
void SKError::calcCovEff(){

  // histogram setup
  hCov = new TH2D("hcov","hcov",Nclass,0,Nclass,Nclass,0,Nclass);
  hCor = new TH2D("hcor","hcor",Nclass,0,Nclass,Nclass,0,Nclass);

  // loops
  for (int i=0; i<Nclass; i++){
    for (int j=0; j<Nclass; j++){
     
      // symmetry
      if (j<i) continue;

      cout<<"finding covariance: "<<i<<","<<j<<endl;
      float cov = arraycov(Efficiency[i],Efficiency[j],Ntoys);
      float cor = arraycor(Efficiency[i],Efficiency[j],Ntoys);
      cout<<"cor: "<<i<<","<<j<<" = "<<cor<<endl;
     
      hCov->SetBinContent(i,j,cov);
      hCor->SetBinContent(i,j,cor);
      //
      if (i!=j){
        hCov->SetBinContent(j,i,cov);
        hCor->SetBinContent(j,i,cor);
      }
    }
  }

  hCor->SetMinimum(-1.);
  hCor->SetMaximum(1.);

  //
  return;
}



///////////////////////////////////////////////////////////////
void SKError::calcCov(){

  // histogram setup
  hCov = new TH2D("hcov","hcov",Nclass,0,Nclass,Nclass,0,Nclass);
  hCor = new TH2D("hcor","hcor",Nclass,0,Nclass,Nclass,0,Nclass);

  // loops
  for (int i=0; i<Nclass; i++){
    for (int j=0; j<Nclass; j++){
     
      // symmetry
      if (j<i) continue;

      cout<<"finding covariance: "<<i<<","<<j<<endl;
      float cov = arraycov(Nevents[i],Nevents[j],Ntoys);
      float cor = arraycor(Nevents[i],Nevents[j],Ntoys);
      cout<<"cor: "<<i<<","<<j<<" = "<<cor<<endl;
     
      hCov->SetBinContent(i,j,cov);
      hCor->SetBinContent(i,j,cor);
      //
      if (i!=j){
        hCov->SetBinContent(j,i,cov);
        hCor->SetBinContent(j,i,cor);
      }
    }
  }

  hCor->SetMinimum(-1.);
  hCor->SetMaximum(1.);

  //
  return;
}



///////////////////////////////////////////////////////////////
void SKError::calcAllEff(int ntoy){

  for (int iclass=0; iclass<Nclass; iclass++){

    // 0th index is nominal
    float eff = calcEff(iclass,ntoy,0); 
    cout<<"eff: "<<iclass<<","<<ntoy<<": "<<eff<<endl;
    Efficiency[iclass][ntoy] = eff;
  }

  return;
}


///////////////////////////////////////////////////////////////
void SKError::addToy(int ntoy){


   // save histos in arrays
   int index = 0;
   //
   for (int ibin=1; ibin<=hEvisNuECCQE->GetNbinsX(); ibin++){
     Nevents[index][ntoy] = hEvisNuECCQE->GetBinContent(ibin);
     index++;
   }
   //
   for (int ibin=1; ibin<=hEvisNuECCOth->GetNbinsX(); ibin++){
     Nevents[index][ntoy] = hEvisNuECCOth->GetBinContent(ibin);
     index++;
   }
   //
   for (int ibin=1; ibin<=hEvisNuMuCCQE->GetNbinsX(); ibin++){
     Nevents[index][ntoy] = hEvisNuMuCCQE->GetBinContent(ibin);
     index++;
   }
   //
   for (int ibin=1; ibin<=hEvisNuMuCCOth->GetNbinsX(); ibin++){
     Nevents[index][ntoy] = hEvisNuMuCCOth->GetBinContent(ibin);
     index++;
   }

   // save histos in TOTAL arrays
   index = 0;
   //
   for (int ibin=1; ibin<=hEvisNuECCQE->GetNbinsX(); ibin++){
     NeventsTotal[index][ntoy] = hEvisNuECCQETot->GetBinContent(ibin);
     index++;
   }
   //
   for (int ibin=1; ibin<=hEvisNuECCOth->GetNbinsX(); ibin++){
     NeventsTotal[index][ntoy] = hEvisNuECCOthTot->GetBinContent(ibin);
     index++;
   }
   //
   for (int ibin=1; ibin<=hEvisNuMuCCQE->GetNbinsX(); ibin++){
     NeventsTotal[index][ntoy] = hEvisNuMuCCQETot->GetBinContent(ibin);
     index++;
   }
   //
   for (int ibin=1; ibin<=hEvisNuMuCCOth->GetNbinsX(); ibin++){
     NeventsTotal[index][ntoy] = hEvisNuMuCCOthTot->GetBinContent(ibin);
     index++;
   }

   // calculate all of the epsilon values!
   calcAllEff(ntoy);

 // 
 return;

}


///////////////////////////////////////////////////////////////
void SKError::resetHistos(){

  //
  hEvisNuECCQE->Reset();
  hEvisNuMuCCQE->Reset();
  hEvisNuECCOth->Reset();
  hEvisNuMuCCOth->Reset();

  //
  return;
}


///////////////////////////////////////////////////////////////
void SKError::addEvent(int nclass, float evis, float weight, bool flgtot){

//  int nclass = getClassMC(nutype, mode, component, evis);

  if (!flgtot){
    //
    if (nclass==1){
      hEvisNuECCQE->Fill(evis,weight);
    }
    if (nclass==2){
      hEvisNuMuCCQE->Fill(evis,weight);
    }
    if (nclass==3){
      hEvisNuECCOth->Fill(evis,weight);
    }
    if (nclass==4){
      hEvisNuMuCCOth->Fill(evis,weight);
    }
  }
  else{
    //
    if (nclass==1){
      hEvisNuECCQETot->Fill(evis,weight);
    }
    if (nclass==2){
      hEvisNuMuCCQETot->Fill(evis,weight);
    }
    if (nclass==3){
      hEvisNuECCOthTot->Fill(evis,weight);
    }
    if (nclass==4){
      hEvisNuMuCCOthTot->Fill(evis,weight);
    }
  }
  //
  return;
}


///////////////////////////////////////////////////////////////
int SKError::getClassMC(int nutype, int mode, int component, float evis){
  
  // core class code: 
  //  1 -> CCQE Nu E
  //  2 -> CCQE Nu Mu
  //  3 -> CCOth Nu E
  //  4 -> CCOth Nu Mu
  //  5 -> NC pi0

//  if ( nsubev==1 && pidlike>=0. && evis>100. && pi0like<0.) return 1;
//  if ( nsubev>1 && pidlike>=0.&& evis>100. && pi0like<0.) return 3;
//  if ( nsubev<=2 && pidlike<=0. && evis>100.) return 2;
//  if ( nsubev>2 && pidlike<=0. && evis>100.) return 4;
//  if ( nsubev==1 && pi0like>0.) return 5;

  // CC and single electron
  if (component==0 && TMath::Abs(mode)<30 ){
    return 1;
  }

  // CC electron and other
  if (component==2 && TMath::Abs(mode)<30){
    return 3;
  }

  // CC muon and single muon
  if (component==1 && TMath::Abs(mode)<30){
    return 2;
  }

  // CC muon and other
  if (component==3 && TMath::Abs(mode)<30){
    return 4;
  }

  // NC
  if (component==4 && TMath::Abs(mode)>=30){
    return 5;
  }

  return 0;
}



///////////////////////////////////////////////////////////////
void SKError::drawAllEff(){

  TH1D* hTmp[NTOYS];

  //
  for (int itoy=0; itoy<Ntoys; itoy++){
    drawSliceEff(itoy);
    hTmp[itoy] = (TH1D*)hSlice->Clone(Form("htmp%d",itoy));
  }
  //
  hTmp[0]->Draw();
  for (int itoy=1; itoy<Ntoys; itoy++){
    hTmp[itoy]->Draw("same");
  }

  //
  return; 
}



///////////////////////////////////////////////////////////////
void SKError::drawAll(){

  TH1D* hTmp[NTOYS];

  //
  for (int itoy=0; itoy<Ntoys; itoy++){
    drawSlice(itoy);
    hTmp[itoy] = (TH1D*)hSlice->Clone(Form("htmp%d",itoy));
  }
  //
  hTmp[0]->Draw();
  for (int itoy=1; itoy<Ntoys; itoy++){
    hTmp[itoy]->Draw("same");
  }

  //
  return; 
}



///////////////////////////////////////////////////////////////
void SKError::drawSlice(int ntoy){

  // reset and fill 
  hSlice->Reset();
  for (int iclass=0; iclass<Nclass; iclass++){
    hSlice->SetBinContent(iclass+1,Nevents[iclass][ntoy]);
    hSlice->SetBinError(iclass+1,0.);
  }

  hSlice->Draw();
  
  //
  return; 
}



///////////////////////////////////////////////////////////////
void SKError::drawSliceEff(int ntoy){

  // reset and fill 
  hSlice->Reset();
  for (int iclass=0; iclass<Nclass; iclass++){
    hSlice->SetBinContent(iclass+1,Efficiency[iclass][ntoy]);
    hSlice->SetBinError(iclass+1,0.);
  }

  hSlice->Draw();
  
  //
  return; 
}



///////////////////////////////////////////////////////////////
void SKError::zeroArrays(){

  // set initial values
  for (int i=0; i<Nclass; i++){
    for (int j=0; j<Nclass; j++){
      Nevents[i][j] = 0.;
      NeventsTotal[i][j] = 0.;
      Efficiency[i][j] = 0.;
    }
  }

  return;
}



///////////////////////////////////////////////////////////////
void SKError::initHistos(){

  // histo binning 
  int    NbinsNuECCQE = 6;
  double BinsNuECCQE[] = {100.,300.,700.,1250.,2000.,5000.,30000.};
  int    NbinsNuECCOth = 3;
  double BinsNuECCOth[] = {100.,1250.,5000.,30000.};
  int    NbinsNuMuCCQE = 7;
  double BinsNuMuCCQE[] = {0.,100.,300.,700.,1250.,2000.,5000.,30000.};
  int    NbinsNuMuCCOth = 3;
  double BinsNuMuCCOth[] = {100.,1250.,5000.,30000.};

  // setup histos
  hEvisNuECCQE = new TH1D("hnueccqe","hnueccqe",NbinsNuECCQE,BinsNuECCQE); 
  hEvisNuECCOth = new TH1D("hnueccoth","hnueccoth",NbinsNuECCOth,BinsNuECCOth); 
  hEvisNuMuCCQE = new TH1D("hnumuccqe","hnumuccqe",NbinsNuMuCCQE,BinsNuMuCCQE); 
  hEvisNuMuCCOth = new TH1D("hnumuccoth","hnumuccoth",NbinsNuMuCCOth,BinsNuMuCCOth); 
  // totals
  hEvisNuECCQETot = new TH1D("hnueccqetot","hnueccqetot",NbinsNuECCQE,BinsNuECCQE); 
  hEvisNuECCOthTot = new TH1D("hnueccothtot","hnueccothtot",NbinsNuECCOth,BinsNuECCOth); 
  hEvisNuMuCCQETot = new TH1D("hnumuccqetot","hnumuccqetot",NbinsNuMuCCQE,BinsNuMuCCQE); 
  hEvisNuMuCCOthTot = new TH1D("hnumuccothtot","hnumuccothtot",NbinsNuMuCCOth,BinsNuMuCCOth); 

  // count all bins
  Nclass = NbinsNuECCQE + NbinsNuECCOth + NbinsNuMuCCQE + NbinsNuMuCCOth;

  // setup slice
  hSlice = new TH1D("hslice","hslice",Nclass,0,Nclass);

  //
  return;
}



///////////////////////////////////////////////////////////////
SKError::SKError(int ntoys){

  // setup the histos
  Ntoys = ntoys;
  initHistos();

}





#endif
