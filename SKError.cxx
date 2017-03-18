#ifndef SKERROR_C
#define SKERROR_C

#include "SKError.h"



///////////////////////////////////////////////////////////////
void SKError::initTN186Errors(){

  // from TN
  double shifts[] = {0.555,1.381,-7.092,1.785,4.6263,11.940
                    ,-7.123,5.891,-10.826
                    ,0.87,1.227,0.094,0.136,-0.037,-0.102,-0.157
                    ,6.406, 2.766, -5.67};
  double fits[] = {1.517,2.045,4.380,2.428,7.912,9.272
                  ,10.187,9.096,11.705
                  ,0.376,0.449,0.582,1.194,0.781,1.334,2.485
                  ,5.771,6.109,11.559};
  int nclasses = 16;
  
  // set values
  for ( int i=0; i<nclasses; i++ ) {
    tn186ShiftError[i]=shifts[i];
    tn186FitError[i]=fits[i];
    tn186TotError[i]=TMath::Sqrt(fits[i]*fits[i] + shifts[i]*shifts[i]);
  }

  // make histos
  hErrorTN186CCQE[0] = (TH1D*)hDiagonalErrorsCCQE[0]->Clone("tn186_nueccqe");
  hErrorTN186CCQE[1] = (TH1D*)hDiagonalErrorsCCQE[1]->Clone("tn186_numccqe");
  hErrorTN186CCOth[0] = (TH1D*)hDiagonalErrorsCCOth[0]->Clone("tn186_nueccoth");
  hErrorTN186CCOth[1] = (TH1D*)hDiagonalErrorsCCOth[1]->Clone("tn186_numccoth");

  // identify bins
  int bins_nue_ccqe[] = {1,2,3,4,5,6};
  int bins_nue_ccoth[] = {7,8,9};
  int bins_num_ccqe[] = {10,11,12,13,14,15,16};
  int bins_num_ccoth[] = {17,18,19};

  // fill in histos

  // ccqe nue
  for (int i=2; i<=hErrorTN186CCQE[0]->GetNbinsX(); i++){
    int index = bins_nue_ccqe[i-2];
    hErrorTN186CCQE[0]->SetBinContent(i,tn186ShiftError[index]);
    hErrorTN186CCQE[0]->SetBinError(i,tn186FitError[index]);
  }
  // ccqe nue
  for (int i=1; i<=hErrorTN186CCQE[1]->GetNbinsX(); i++){
    int index = bins_num_ccqe[i-1];
    hErrorTN186CCQE[1]->SetBinContent(i,tn186ShiftError[index]);
    hErrorTN186CCQE[1]->SetBinError(i,tn186FitError[index]);
  }
  for (int i=1; i<=hErrorTN186CCOth[0]->GetNbinsX(); i++){
    int index = bins_nue_ccoth[i-1];
    hErrorTN186CCOth[0]->SetBinContent(i,tn186ShiftError[index]);
    hErrorTN186CCOth[0]->SetBinError(i,tn186FitError[index]);
  }

  //
  return;
}



///////////////////////////////////////////////////////////////
double SKError::getMaxError(TH1D* hh){
  double errmax = 0.;
  for ( int i=1; i<=hh->GetNbinsX(); i++ ) {
    double err = hh->GetBinError(i);
    if (err>errmax){
      errmax = err;
    }
  }
  cout<<"max err: "<<errmax<<endl;
  return errmax;
}



///////////////////////////////////////////////////////////////
std::vector<TLatex*> SKError::getBinLabels(TH1D*hh){

  std::vector<TLatex*> vlabel;
  double xstart = -0.5;
  for (int ih=0; ih<4; ih++){
    for (int i=1; i<=hh->GetNbinsX(); i++){
      double xpos = xstart + i;
      double offset = -0.5;
      double evismin = hh->GetXaxis()->GetBinLowEdge(i)/1000.;
      double evismax = hh->GetXaxis()->GetBinWidth(i)/1000. + evismin;
      vlabel.push_back(new TLatex(offset,xpos,Form("[%1.1f - %1.1f]",evismin,evismax)));
      vlabel.at(i-1)->SetTextAlign(32);
      vlabel.at(i-1)->SetTextSize(0.08);
    }
  }

  return vlabel;
}



///////////////////////////////////////////////////////////////
void SKError::drawDiagonals(){

  // calc first
  calcDiagonals();

  // canvas setup
  TCanvas* cc = new TCanvas("cc","cc",800,600);
  cc->Divide(1,2);

  // histogram properties
  hDiagonalErrorsCCQE[0]->SetLineColor(kBlue);
  hDiagonalErrorsCCQE[0]->SetFillColor(kBlue);
  hDiagonalErrorsCCQE[0]->SetLineWidth(3);
  hDiagonalErrorsCCQE[0]->SetMarkerStyle(0);
  hDiagonalErrorsCCQE[1]->SetLineColor(kRed);
  hDiagonalErrorsCCQE[1]->SetFillColor(kRed);
  hDiagonalErrorsCCQE[1]->SetLineWidth(3);
  hDiagonalErrorsCCQE[1]->SetMarkerStyle(0);
  hDiagonalErrorsCCOth[0]->SetLineColor(kBlue);
  hDiagonalErrorsCCOth[0]->SetFillColor(kBlue);
  hDiagonalErrorsCCOth[0]->SetLineWidth(3);
  hDiagonalErrorsCCOth[0]->SetMarkerStyle(0);
  hDiagonalErrorsCCOth[1]->SetLineColor(kRed);
  hDiagonalErrorsCCOth[1]->SetFillColor(kRed);
  hDiagonalErrorsCCOth[1]->SetLineWidth(3);
  hDiagonalErrorsCCOth[1]->SetMarkerStyle(0);

  // find ranges
  double margin = 1.1;
  double range_nue_ccqe = TMath::Max( TMath::Abs(hDiagonalErrorsCCQE[0]->GetMaximum()),
                                      TMath::Abs(hDiagonalErrorsCCQE[0]->GetMinimum()));
  double range_numu_ccqe = TMath::Max( TMath::Abs(hDiagonalErrorsCCQE[1]->GetMaximum()),
                                      TMath::Abs(hDiagonalErrorsCCQE[1]->GetMinimum()));
  double range_ccqe = margin*TMath::Max( range_nue_ccqe, range_numu_ccqe);
  range_ccqe += getMaxError(hDiagonalErrorsCCQE[0]);

  cout<<"range: "<<range_ccqe<<endl; 
  //
  double range_nue_ccoth = TMath::Max( TMath::Abs(hDiagonalErrorsCCOth[0]->GetMaximum()),
                                      TMath::Abs(hDiagonalErrorsCCOth[0]->GetMinimum()));
  double range_numu_ccoth = TMath::Max( TMath::Abs(hDiagonalErrorsCCOth[1]->GetMaximum()),
                                      TMath::Abs(hDiagonalErrorsCCOth[1]->GetMinimum()));
  double range_ccoth = margin*TMath::Max( range_nue_ccoth, range_numu_ccoth);
  range_ccoth += getMaxError(hDiagonalErrorsCCOth[0]);

  cout<<"range: "<<range_ccoth<<endl; 
  //
  hDiagonalErrorsCCQE[0]->SetMinimum(-1.*range_ccqe);
  hDiagonalErrorsCCQE[0]->SetMaximum(1.*range_ccqe);
  //
  hDiagonalErrorsCCQE[1]->SetMinimum(-1.*range_ccqe);
  hDiagonalErrorsCCQE[1]->SetMaximum(1.*range_ccqe);
  //
  hDiagonalErrorsCCOth[0]->SetMinimum(-1.*range_ccoth);
  hDiagonalErrorsCCOth[0]->SetMaximum(1.*range_ccoth);
  //
  hDiagonalErrorsCCOth[1]->SetMinimum(-1.*range_ccoth);
  hDiagonalErrorsCCOth[1]->SetMaximum(1.*range_ccoth);


  // get bin labels
  vector<TLatex*> vlabel_ccqe = getBinLabels(hEvisNuMuCCQE);
  vector<TLatex*> vlabel_ccoth = getBinLabels(hEvisNuECCOth);

  // draw
  cc->cd(1);
  for (int ibin=1; ibin<=hDiagonalErrorsCCQE[0]->GetNbinsX(); ibin++){
    hDiagonalErrorsCCQE[0]->GetXaxis()->SetBinLabel(ibin,vlabel_ccqe.at(ibin-1)->GetTitle());
  }
  hDiagonalErrorsCCQE[0]->GetXaxis()->SetLabelSize(0.08);
  hDiagonalErrorsCCQE[0]->Draw("E");
  hDiagonalErrorsCCQE[1]->Draw("sameE");
  //
  cc->cd(2);
  for (int ibin=1; ibin<=hDiagonalErrorsCCOth[0]->GetNbinsX(); ibin++){
    hDiagonalErrorsCCOth[0]->GetXaxis()->SetBinLabel(ibin,vlabel_ccqe.at(ibin-1)->GetTitle());
  }
  hDiagonalErrorsCCOth[0]->GetXaxis()->SetLabelSize(0.08);
  hDiagonalErrorsCCOth[0]->Draw("E");
  hDiagonalErrorsCCOth[1]->Draw("sameE");

  //
  return;

}


///////////////////////////////////////////////////////////////
void SKError::calcDiagonals(){

  // identify bins
  int bins_nue_ccqe[] = {1,2,3,4,5,6};
  int bins_nue_ccoth[] = {7,8,9};
  int bins_num_ccqe[] = {10,11,12,13,14,15,16};
  int bins_num_ccoth[] = {17,18,19};

  // ccqe nue
  for (int i=hDiagonalErrorsCCQE[0]->GetNbinsX(); i>1; i--){
    int nbin = bins_nue_ccqe[i-2];
    hDiagonalErrorsCCQE[0]->SetBinContent(i,hDiagonalErrors->GetBinContent(nbin));
    hDiagonalErrorsCCQE[0]->SetBinError(i,hDiagonalErrors->GetBinError(nbin));
  }

  // ccqe numu
  for (int i=hDiagonalErrorsCCQE[1]->GetNbinsX(); i>0; i--){
    // fill numu ccqe bins
    int nbin = bins_num_ccqe[i-1];
    hDiagonalErrorsCCQE[1]->SetBinContent(i,hDiagonalErrors->GetBinContent(nbin));
    hDiagonalErrorsCCQE[1]->SetBinError(i,hDiagonalErrors->GetBinError(nbin));
  }

  // ccother nue
  for (int i=hDiagonalErrorsCCOth[0]->GetNbinsX(); i>0; i--){
    // fill numu ccqe bins
    int nbin = bins_nue_ccoth[i-1];
    hDiagonalErrorsCCOth[0]->SetBinContent(i,hDiagonalErrors->GetBinContent(nbin));
    hDiagonalErrorsCCOth[0]->SetBinError(i,hDiagonalErrors->GetBinError(nbin));
  }

  // ccother numu
  for (int i=hDiagonalErrorsCCOth[1]->GetNbinsX(); i>0; i--){
    // fill numu ccqe bins
    int nbin = bins_num_ccoth[i-1];
    hDiagonalErrorsCCOth[1]->SetBinContent(i,hDiagonalErrors->GetBinContent(nbin));
    hDiagonalErrorsCCOth[1]->SetBinError(i,hDiagonalErrors->GetBinError(nbin));
  }

  return;
}



///////////////////////////////////////////////////////////////
double SKError::calcFitError(int iclass){
    cout<<"Fit error for class: "<<iclass<<endl;
    return TMath::Sqrt(arrayvarD(DelEfficiency[iclass],Ntoys));
}



///////////////////////////////////////////////////////////////
void SKError::calcErrors(){
    for (int iclass=0; iclass<Nclass; iclass++){
      DelEffShiftError[iclass] = calcShiftError(iclass);
      absDelEffShiftError[iclass] = TMath::Abs(calcShiftError(iclass));
      DelEffFitError[iclass] = calcFitError(iclass);
    }
    vShiftErrors = new TVectorD(Nclass,DelEffShiftError);
    return;
}



///////////////////////////////////////////////////////////////
void SKError::printErrors(){
    for (int iclass=0; iclass<Nclass; iclass++){
      cout<<"--- Class "<<iclass<<"  ---"<<endl;
      cout<<"  Fit: "<<DelEffFitError[iclass]*100.<<"%";
      cout<<"  Shift: "<<DelEffShiftError[iclass]*100.<<"%";
      double toterr = TMath::Sqrt(DelEffFitError[iclass]*DelEffFitError[iclass]
                                 +DelEffShiftError[iclass]*DelEffShiftError[iclass]);
      cout<<"  Total: "<<toterr*100.<<"%"<<endl;
    }
}



///////////////////////////////////////////////////////////////
void SKError::printEffDist(const char* plotdir){
     TCanvas* cc = new TCanvas("cc","cc",700,600);
     cc->GetPad(0)->SetLeftMargin(0.12);
     cc->GetPad(0)->SetRightMargin(0.12);
     cc->GetPad(0)->SetBottomMargin(0.15);
     TString pdir = plotdir;
     for (int i=0; i<Nclass; i++){
       TString pname = pdir.Data();
       pname.Append(Form("throw_dist_class_%d.png",i));
       drawEffDist(i);
       cc->Print(pname.Data());
     }
     return;
}



///////////////////////////////////////////////////////////////
void SKError::saveErrors(const char* filename){
    TFile* fout = new TFile(filename,"RECREATE");
    hCov->Write();
    hCor->Write();
    vShiftErrors->Write("vshift");
    cout<<"writing: "<<filename<<endl;
    fout->Close();
    return;
}



///////////////////////////////////////////////////////////////
double SKError::calcShiftError(int iclass){
    cout<<"Shift error for class: "<<iclass<<endl;
    return arraymeanD(DelEfficiency[iclass],Ntoys);
}



///////////////////////////////////////////////////////////////
void SKError::drawVertLines(){
  
  // line properties have been set in initHistos()
  
  int nlines = 3;

  for ( int i=0; i<nlines; i++ ) {
    lineVert[i]->Draw("same");
  }

}



///////////////////////////////////////////////////////////////
void SKError::drawHorizLines(){
  
  // line properties have been set in initHistos()
  
  int nlines = 3;

  for ( int i=0; i<nlines; i++ ) {
    lineHorz[i]->Draw("same");
  }

}



///////////////////////////////////////////////////////////////
void SKError::makeBinLabels(){

  // bin labels
  TH1D* histos[]={hEvisNuECCQE,hEvisNuECCOth,hEvisNuMuCCQE,hEvisNuMuCCOthTot};
  double xstart = -0.5;
  int iclass = 0;
  for (int ih=0; ih<4; ih++){
    for (int i=1; i<=histos[ih]->GetNbinsX(); i++){
      double xpos = xstart + i;
      double offset = -0.5;
      double evismin = histos[ih]->GetXaxis()->GetBinLowEdge(i)/1000.;
      double evismax = histos[ih]->GetXaxis()->GetBinWidth(i)/1000. + evismin;
      labelHorz[iclass] = new TLatex(xpos,offset,Form("[%1.1f - %1.1f]",evismin,evismax));
      labelVert[iclass] = new TLatex(offset,xpos,Form("[%1.1f - %1.1f]",evismin,evismax));
      labelHorz[iclass]->SetTextAngle(90);
      labelHorz[iclass]->SetTextAlign(32);
      labelVert[iclass]->SetTextAlign(32);
      labelHorz[iclass]->SetTextSize(0.020);
      labelVert[iclass]->SetTextSize(0.020);
      iclass++;
    }
    xstart += (double)histos[ih]->GetNbinsX();
  }

  // sector labels
  TString sector[] = {"CCQE","CCOth.","CCQE","CCOth."};
  TString nulabel[] = {"#nu_{e}","#nu_{e}","#nu_{#mu}","#nu_{#mu}"};
  double labelpos[]= {3.,7.5,12.5,17.5};
  double offset = -3.5;
  double nuoffset = -4.5;
  for (int i=0; i<4; i++){
   sectorLabelHorz[i] = new TLatex(labelpos[i],offset,sector[i].Data()); 
   nuLabelHorz[i] = new TLatex(labelpos[i],nuoffset,nulabel[i].Data());
   sectorLabelVert[i] = new TLatex(offset,labelpos[i],sector[i].Data()); 
   nuLabelVert[i] = new TLatex(nuoffset,labelpos[i],nulabel[i].Data());
   sectorLabelVert[i]->SetTextAngle(90); 
   sectorLabelHorz[i]->SetTextAlign(22); 
   sectorLabelVert[i]->SetTextAlign(22); 
   nuLabelHorz[i]->SetTextAlign(22); 
   nuLabelVert[i]->SetTextAlign(22); 
   sectorLabelVert[i]->SetTextSize(0.03); 
   sectorLabelHorz[i]->SetTextSize(0.03); 
  }

  return;

}



///////////////////////////////////////////////////////////////
void SKError::drawBinLabels(){

  // bin labels
  TH1D* histos[]={hEvisNuECCQE,hEvisNuECCOth,hEvisNuMuCCQE,hEvisNuMuCCOthTot};
  int iclass = 0;
  for (int ih=0; ih<4; ih++){
    for (int i=1; i<=histos[ih]->GetNbinsX(); i++){
      labelHorz[iclass]->Draw("same");
      labelVert[iclass]->Draw("same");
      iclass++;
    }
  }

  // sector labels
  for (int i=0; i<4; i++){
   sectorLabelVert[i]->Draw("same"); 
   sectorLabelHorz[i]->Draw("same"); 
   nuLabelHorz[i]->Draw("same");
   nuLabelVert[i]->Draw("same");
  }

  return;
}



///////////////////////////////////////////////////////////////
void SKError::drawCor(){
  
  TCanvas *cc = new TCanvas("cc","cc",750,700);
  cc->GetPad(0)->SetLeftMargin(0.22);
  cc->GetPad(0)->SetRightMargin(0.15);
  cc->GetPad(0)->SetTopMargin(0.1);
  cc->GetPad(0)->SetBottomMargin(0.22);

  hCor->Draw("colz");
  hCor->GetXaxis()->SetNdivisions(0);
  hCor->GetYaxis()->SetNdivisions(0);
  hCor->GetZaxis()->SetTitle("Correlation");
  hCor->GetZaxis()->SetTitleSize(0.05);
  hCor->GetZaxis()->SetTitleOffset(0.8);
  hCor->GetZaxis()->CenterTitle(1);


  drawVertLines();
  drawHorizLines();
  drawBinLabels();


  return;
}



///////////////////////////////////////////////////////////////
void SKError::drawCov(){

  hCov->Draw("colz");
  TCanvas *cc = new TCanvas("cc","cc",750,700);
  cc->GetPad(0)->SetLeftMargin(0.12);
  cc->GetPad(0)->SetRightMargin(0.12);
  cc->GetPad(0)->SetTopMargin(0.12);
  cc->GetPad(0)->SetBottomMargin(0.12);
  for (int i=0; i<NLINES; i++){
    lineHorz[i]->Draw("same");
    lineVert[i]->Draw("same");
  }
  return;
}




///////////////////////////////////////////////////////////////
// calcutlate efficiency
// here this is the ratio of modified to original eff
double SKError::calcEff(int nclass, int ntoy){

  // don't divide by zero
  if (NeventsTotal[nclass][0]==0) return 0.;

  // get this efficiency  
  double eff = Nevents[nclass][ntoy] / NeventsTotal[nclass][ntoy];
  cout<<"eff: "<<eff<<endl;
  double eff0 = Nevents[nclass][0] / NeventsTotal[nclass][0];
  cout<<"eff0: "<<eff0<<endl;
  
  if (eff0 == 0) return 0.;

  return ((eff/eff0) - 1.); 

}



///////////////////////////////////////////////////////////////
// calcutlate efficiency
// here the modified eff. is calc w.r.t. modified total
double SKError::calcDelEff(int nclass, int ntoy){

  // don't divide by zero
  if (NeventsTotal[nclass][0]==0) return 0.;
  if (NeventsTotal[nclass][ntoy]==0) return 0.;

  // get this efficiency  
  double eff = Nevents[nclass][ntoy] / NeventsTotal[nclass][ntoy];
  cout<<"eff: "<<eff<<endl;
  double eff0 = Nevents[nclass][0] / NeventsTotal[nclass][0];
  cout<<"eff0: "<<eff0<<endl;
  
  if (eff0 == 0) return 0.;

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

  gScat->Draw("a*");

  return;
}



///////////////////////////////////////////////////////////////
void SKError::calcCov(int vartype){


  // loops
  for (int i=0; i<Nclass; i++){

    // fill diagonal errors
    
    // efficiency errors
    if ( vartype==0 ) {
      hDiagonalErrors->SetBinContent(i+1,arraymeanD(DelEfficiency[i],Ntoys));
      hDiagonalErrors->SetBinError(i+1,TMath::Sqrt(arrayvarD(DelEfficiency[i],Ntoys)));
      double cov = arraycovD(DelEfficiency[i],DelEfficiency[i],Ntoys);
      double cor = arraycorD(DelEfficiency[i],DelEfficiency[i],Ntoys);
      hCov->SetBinContent(i+1,i+1,cov);
      hCor->SetBinContent(i+1,i+1,cor);
    }

    // NSK errors
    else if ( vartype==1 ){
      hDiagonalErrors->SetBinContent(i+1,arraymeanD(Nevents[i],Ntoys));
      hDiagonalErrors->SetBinError(i+1,TMath::Sqrt(arrayvarD(Nevents[i],Ntoys)));       
      double cov = arraycovD(Nevents[i],Nevents[i],Ntoys);
      double cor = arraycorD(Nevents[i],Nevents[i],Ntoys);
      hCov->SetBinContent(i+1,i+1,cov);
      hCor->SetBinContent(i+1,i+1,cor);
    }

    // fill off-diagonal errors
    for (int j=0; j<Nclass; j++){
     
      // symmetry
      if (j<=i) continue;

      // off-diagonal elements
      if ( vartype==0 ) {
        double cov = arraycovD(DelEfficiency[i],DelEfficiency[j],Ntoys);
        double cor = arraycorD(DelEfficiency[i],DelEfficiency[j],Ntoys);
        hCov->SetBinContent(i+1,j+1,cov);
        hCor->SetBinContent(i+1,j+1,cor);
        hCov->SetBinContent(j+1,i+1,cov);
        hCor->SetBinContent(j+1,i+1,cor);
      }
      else if ( vartype==1 ){
        double cov = arraycovD(Nevents[i],Nevents[j],Ntoys);
        double cor = arraycorD(Nevents[i],Nevents[j],Ntoys);
        hCov->SetBinContent(i+1,j+1,cov);
        hCor->SetBinContent(i+1,j+1,cor);
        hCov->SetBinContent(j+1,i+1,cov);
        hCor->SetBinContent(j+1,i+1,cor);
      }
    }
  }

  //scales
  hCor->SetMinimum(-1.11);
  hCor->SetMaximum(1.11);

  // draw correlations
  drawCor();

  return;
}



///////////////////////////////////////////////////////////////
void SKError::calcAllDelEff(int ntoy, int effdef){

  // find efficiency for each class
  for (int iclass=0; iclass<Nclass; iclass++){

    double eff = 0.;
    if (effdef==0){
      eff = calcDelEff(iclass,ntoy); 
      DelEfficiency[iclass][ntoy] = eff;
    }

    // use nominal value (instead of alpha-modified, see TN-186)
    // in the denominator:
    else if (effdef==1){
      eff = calcEff(iclass,ntoy); 
      DelEfficiency[iclass][ntoy] = eff;
    }

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
   //;
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
   calcAllDelEff(ntoy);

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
  hEvisNuECCQETot->Reset();
  hEvisNuMuCCQETot->Reset();
  hEvisNuECCOthTot->Reset();
  hEvisNuMuCCOthTot->Reset();

  return;
}


///////////////////////////////////////////////////////////////
void SKError::addEvent(int nclass, double evis, double weight, bool flgtot){


  // passes core
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

  // total
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
int SKError::getClassMC(int nutype, int mode, int component,
                        double evis, int nsubev, double towall, double wall){
  
  // class code: 
  //
  //  0 -> Undefined
  //  1 -> CCQE Nu E
  //  2 -> CCQE Nu Mu
  //  3 -> CCOth Nu E
  //  4 -> CCOth Nu Mu
  //  5 -> NC pi0

  double towall_nuecut = 170.;
  double towall_numucut = 250.;
  double wall_nuecut = 80.;
  double wall_numucut = 50.;

  
  // CC and single electron
  if (component==0 && TMath::Abs(mode)<30 && evis>100. && nutype==12
      && nsubev==1 && towall>towall_nuecut && wall>wall_nuecut){
    return 1;
  }

  // CC electron and other
  if (component==2 && TMath::Abs(mode)<30 && evis>100. && nutype==12
      && nsubev>1 && towall>towall_nuecut && wall>wall_nuecut ){
    return 3;
  }

  // CC and single muon
  if (component==1 && TMath::Abs(mode)<30 && nsubev==2 && nutype==14
      && evis>30. && towall>towall_numucut && wall>wall_numucut  ){
    return 2;
  }

  // CC muon and other
  if (component==3 && TMath::Abs(mode)<30 && nsubev>2 && evis>30.
       && towall>towall_numucut && wall>wall_numucut && nutype==14){
    return 4;
  }

  // NC
  if (component==4 && TMath::Abs(mode)>=30 && evis>100. && nsubev==1
       && towall>towall_nuecut && wall>wall_nuecut ){
    return 5;
  }
  

  /*
  // CC and single electron
  if (component==0 && TMath::Abs(mode)<30 && evis>100. && nutype==12
      && nsubev==1 && towall>towall_nuecut && wall>wall_nuecut){
    return 1;
  }

  // CC electron and other
  if (component==2 && TMath::Abs(mode)<30 && evis>100. && nutype==12
      && nsubev==1 && towall>towall_nuecut && wall>wall_nuecut ){
    return 3;
  }

  // CC and single muon
  if (component==1 && TMath::Abs(mode)<30 && nsubev==2 && nutype==14
      && evis>30. && towall>towall_numucut && wall>wall_numucut  ){
    return 2;
  }

  // CC and other
  if (component==3 && TMath::Abs(mode)<30 && nsubev==2 && evis>30.
       && towall>towall_numucut && wall>wall_numucut && nutype==14){
    return 4;
  }

  // NC
  if (component==4 && TMath::Abs(mode)>=30 && evis>100. && nsubev==1
       && towall>towall_nuecut && wall>wall_nuecut ){
    return 5;
  }
  */

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
void SKError::drawDist(int iclass){

   // setup histogram range
   double xmin = arrayminD(Nevents[iclass],Ntoys);
   cout<<"xmin: "<<xmin<<endl;
   double xmax = arraymaxD(Nevents[iclass],Ntoys);
   cout<<"xmax: "<<xmax<<endl;
   double width = (xmax-xmin);
   xmin = xmin-(width/2.);
   xmax = xmax+(width/2.);

   hdist=new TH1D("hdist","hdist",30,xmin,xmax);

   for (int i=0; i<Ntoys; i++){
     hdist->Fill(Nevents[iclass][i]);
   }

   hdist->SetLineWidth(3);
   hdist->GetXaxis()->SetTitle("#Delta #epsilon");
   hdist->GetXaxis()->SetNdivisions(5);
   hdist->Draw();

   return;
}



///////////////////////////////////////////////////////////////
void SKError::drawEffDist(int iclass){
 
   // setup histogram range
   double xmin = arrayminD(DelEfficiency[iclass],Ntoys);
   cout<<"xmin: "<<xmin<<endl;
   double xmax = arraymaxD(DelEfficiency[iclass],Ntoys);
   cout<<"xmax: "<<xmax<<endl;
   double width = (xmax-xmin);
   xmin = xmin-(width/2.);
   xmax = xmax+(width/2.);

   hdist=new TH1D("hdist","hdist",30,xmin,xmax);

   for (int i=0; i<Ntoys; i++){
     hdist->Fill(DelEfficiency[iclass][i]);
   }

   hdist->SetLineWidth(3);
   hdist->SetLineColor(kBlue);
   hdist->GetXaxis()->SetTitle("#Delta #varepsilon");
   hdist->GetYaxis()->SetTitle("# of throws");
   hdist->GetYaxis()->SetTitleSize(0.06);
   hdist->GetXaxis()->SetTitleSize(0.06);
   hdist->GetXaxis()->SetTitleOffset(0.8);
   hdist->GetYaxis()->SetTitleOffset(0.8);
   hdist->GetXaxis()->SetNdivisions(5);
   hdist->Draw();

   // lines
   TLine* lmean = new TLine(hdist->GetMean(),0,hdist->GetMean(),hdist->GetMaximum());
   lmean->SetLineWidth(3);
   lmean->SetLineColor(kRed);
   lmean->SetLineStyle(2);
   lmean->Draw("same");

   TLine* lzero = new TLine(0,0,0,hdist->GetMaximum());
   lzero->SetLineWidth(2);
   lzero->Draw("same");

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
void SKError::drawSliceTot(int ntoy){

  // reset and fill 
  hSlice->Reset();
  for (int iclass=0; iclass<Nclass; iclass++){
    hSlice->SetBinContent(iclass+1,NeventsTotal[iclass][ntoy]);
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
    hSlice->SetBinContent(iclass+1,DelEfficiency[iclass][ntoy]);
    hSlice->SetBinError(iclass+1,0.);
  }
  hSlice->SetMinimum(-0.3);
  hSlice->SetMaximum(0.3);

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
      DelEfficiency[i][j] = 0.;
    }
  }

  return;
}



///////////////////////////////////////////////////////////////
void SKError::initHistos(int ibinning){

  // histo binning (patrick :) )
  if (ibinning==1){
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
  }

  if (ibinning==0){
  // histo binning (andy :) )
  int    NbinsNuECCQE = 1;
  double BinsNuECCQE[] = {100,30000.};
  int    NbinsNuECCOth = 1;
  double BinsNuECCOth[] = {100.,30000.};
  int    NbinsNuMuCCQE = 3;
  double BinsNuMuCCQE[] = {0.,700.,5000.,30000.};
  int    NbinsNuMuCCOth = 1;
  double BinsNuMuCCOth[] = {0,30000.};
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
  }



  if (ibinning==2){
  // alt. binnings
  
  int nbins = 10;
  double Evismax = 2000;
  hEvisNuECCQE = new TH1D("hnueccqe","hnueccqe",nbins,0,Evismax); 
  hEvisNuECCOth = new TH1D("hnueccoth","hnueccoth",nbins,0,Evismax); 
  hEvisNuMuCCQE = new TH1D("hnumuccqe","hnumuccqe",nbins,0,Evismax); 
  hEvisNuMuCCOth = new TH1D("hnumuccoth","hnumuccoth",nbins,0,Evismax); 

  // totals
  hEvisNuECCQETot = new TH1D("hnueccqetot","hnueccqetot",nbins,0,Evismax); 
  hEvisNuECCOthTot = new TH1D("hnueccothtot","hnueccothtot",nbins,0,Evismax);
  hEvisNuMuCCQETot = new TH1D("hnumuccqetot","hnumuccqetot",nbins,0,Evismax);
  hEvisNuMuCCOthTot = new TH1D("hnumuccothtot","hnumuccothtot",nbins,0,Evismax);
  }


  // count all bins
  Nclass = hEvisNuECCQE->GetNbinsX() + hEvisNuECCOth->GetNbinsX()
           +  hEvisNuMuCCQE->GetNbinsX() + hEvisNuMuCCOth->GetNbinsX();

  // setup slice
  hSlice = new TH1D("hslice","hslice",Nclass,0,Nclass);

  // setup lines for covariance matricies
  int nlines = 3;
  TH1D* hclasses[] = {hEvisNuECCQE,hEvisNuECCOth,hEvisNuMuCCQE};
  double minval = -3.;
  double maxval = Nclass;
  for (int i=0; i<nlines; i++){
    int nbinstot = hclasses[i]->GetNbinsX();
    lineVal[i] = (double)nbinstot;
    if (i>0){
      lineVal[i]+=lineVal[i-1];
    }
    lineHorz[i] = new TLine(minval,lineVal[i],maxval,lineVal[i]);
    lineVert[i] = new TLine(lineVal[i],minval,lineVal[i],maxval);
    lineVert[i]->SetLineWidth(3);
    lineHorz[i]->SetLineWidth(3);
    if (i!=1){
      lineVert[i]->SetLineStyle(2);
      lineHorz[i]->SetLineStyle(2);
    }
  }


  // output histogram setup
  hCov = new TH2D("hcov","hcov",Nclass,0,Nclass,Nclass,0,Nclass);
  hCov = new TH2D("hcov","hcov",Nclass,0,Nclass,Nclass,0,Nclass);
  hCor = new TH2D("hcor","hcor",Nclass,0,Nclass,Nclass,0,Nclass);
  hDiagonalErrors = new TH1D("hdiag","hdiag",Nclass,0,Nclass);
  int nbins_ccqe = TMath::Max(hEvisNuECCQE->GetNbinsX(),hEvisNuMuCCQE->GetNbinsX());
  int nbins_ccoth = TMath::Max(hEvisNuECCOth->GetNbinsX(),hEvisNuMuCCOth->GetNbinsX());
  for ( int i=0; i<2; i++ ) {
    hDiagonalErrorsCCQE[i] = new TH1D("hdiag_ccqe","hdiag_ccqe",nbins_ccqe,0,nbins_ccqe);
    hDiagonalErrorsCCOth[i] = new TH1D("hdiag_ccoth","hdiag_ccoth",nbins_ccoth,0,nbins_ccoth);
  }
  hCor = new TH2D("hcor","hcor",Nclass,0,Nclass,Nclass,0,Nclass);



  // make labels for bins, sections
  makeBinLabels();

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
