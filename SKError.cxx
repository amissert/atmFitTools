#ifndef SKERROR_C
#define SKERROR_C

#include "SKError.h"

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



  // bin labels
  TH1D* histos[]={hEvisNuECCQE,hEvisNuECCOth,hEvisNuMuCCQE,hEvisNuMuCCOthTot};
  double xstart = -0.5;
  for (int ih=0; ih<4; ih++){
    for (int i=1; i<=histos[ih]->GetNbinsX(); i++){

      double xpos = xstart + i;
      double offset = -0.5;
      double evismin = histos[ih]->GetXaxis()->GetBinLowEdge(i)/1000.;
      double evismax = histos[ih]->GetXaxis()->GetBinWidth(i)/1000. + evismin;

      labelHorz[i] = new TLatex(xpos,offset,Form("[%1.1f - %1.1f]",evismin,evismax));
      labelVert[i] = new TLatex(offset,xpos,Form("[%1.1f - %1.1f]",evismin,evismax));
      labelHorz[i]->SetTextAngle(90);
      labelHorz[i]->SetTextAlign(32);
      labelVert[i]->SetTextAlign(32);
      labelHorz[i]->SetTextSize(0.020);
      labelVert[i]->SetTextSize(0.020);
      labelHorz[i]->Draw("same");
      labelVert[i]->Draw("same");
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
   sectorLabelVert[i]->Draw("same"); 
   sectorLabelHorz[i]->Draw("same"); 
   nuLabelHorz[i]->Draw("same");
   nuLabelVert[i]->Draw("same");

  }


  for (int i=0; i<NLINES; i++){
    lineHorz[i]->Draw("same");
    lineVert[i]->Draw("same");
  }
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
void SKError::calcCovDelEff(){

  // histogram setup
  hCov = new TH2D("hcov","hcov",Nclass,0,Nclass,Nclass,0,Nclass);
  hCor = new TH2D("hcor","hcor",Nclass,0,Nclass,Nclass,0,Nclass);
  hDiagonalErrors = new TH1D("hdiag","hdiag",Nclass,0,Nclass);

  // loops
  for (int i=0; i<Nclass; i++){
    for (int j=0; j<Nclass; j++){
     
      // symmetry
      if (j<i) continue;
      cout<<"finding covariance: "<<i<<","<<j<<endl;
      double cov = arraycovD(DelEfficiency[i],DelEfficiency[j],Ntoys);
      double cor = arraycorD(DelEfficiency[i],DelEfficiency[j],Ntoys);
      cout<<"cor: "<<i<<","<<j<<" = "<<cor<<endl;
      hCov->SetBinContent(i+1,j+1,cov);
      hCor->SetBinContent(i+1,j+1,cor);

      //
      if (i!=j){
        hCov->SetBinContent(j+1,i+1,cov);
        hCor->SetBinContent(j+1,i+1,cor);
      }

      // fill diagonal errors
      if (i==j){
        hDiagonalErrors->SetBinContent(i+1,arraymeanD(DelEfficiency[i],Ntoys));
        hDiagonalErrors->SetBinError(i+1,TMath::Sqrt(arrayvarD(DelEfficiency[i],Ntoys)));
      }

    }
  }

  hCor->SetMinimum(-1.11);
  hCor->SetMaximum(1.11);


//  hCor->Draw("colz");
  drawCor();

  //
  return;
}


///////////////////////////////////////////////////////////////
void SKError::calcCovEff(){

  // histogram setup
  hCov = new TH2D("hcov","hcov",Nclass,0,Nclass,Nclass,0,Nclass);
 
  // histogram setup
  hCov = new TH2D("hcov","hcov",Nclass,0,Nclass,Nclass,0,Nclass);
  hCor = new TH2D("hcor","hcor",Nclass,0,Nclass,Nclass,0,Nclass);
  hDiagonalErrors = new TH1D("hdiag","hdiag",Nclass,0,Nclass); hCor = new TH2D("hcor","hcor",Nclass,0,Nclass,Nclass,0,Nclass);

  // loops
  for (int i=0; i<Nclass; i++){
    for (int j=0; j<Nclass; j++){
     
      // symmetry
      if (j<i) continue;
      cout<<"finding covariance: "<<i<<","<<j<<endl;
      double cov = arraycovD(DelEfficiency[i],DelEfficiency[j],Ntoys);
      double cor = arraycorD(DelEfficiency[i],DelEfficiency[j],Ntoys);
      cout<<"cor: "<<i<<","<<j<<" = "<<cor<<endl;
      hCov->SetBinContent(i+1,j+1,cov);
      hCor->SetBinContent(i+1,j+1,cor);

      //
      if (i!=j){
        hCov->SetBinContent(j+1,i+1,cov);
        hCor->SetBinContent(j+1,i+1,cor);
      }

      // fill diagonal errors
      if (i==j){
        hDiagonalErrors->SetBinContent(i+1,arraymeanD(DelEfficiency[i],Ntoys));
        hDiagonalErrors->SetBinError(i+1,TMath::Sqrt(arrayvarD(DelEfficiency[i],Ntoys)));
      }

//      else{ 
//        hCov->SetBinContent(j+1,i+1,cov);
//        hCor->SetBinContent(j+1,i+1,0);
//      }

    }
  }

  hCor->SetMinimum(-1.11);
  hCor->SetMaximum(1.11);


  drawCor();
  //
  return;
}



///////////////////////////////////////////////////////////////
void SKError::calcCov(){

  // histogram setup
  hCov = new TH2D("hcov","hcov",Nclass,0,Nclass,Nclass,0,Nclass);
 
  // histogram setup
  hCov = new TH2D("hcov","hcov",Nclass,0,Nclass,Nclass,0,Nclass);
  hCor = new TH2D("hcor","hcor",Nclass,0,Nclass,Nclass,0,Nclass);
  hDiagonalErrors = new TH1D("hdiag","hdiag",Nclass,0,Nclass); hCor = new TH2D("hcor","hcor",Nclass,0,Nclass,Nclass,0,Nclass);

  // loops
  for (int i=0; i<Nclass; i++){
    for (int j=0; j<Nclass; j++){
     
      // symmetry
      if (j<i) continue;

      cout<<"finding covariance: "<<i<<","<<j<<endl;
      double cov = arraycovD(Nevents[i],Nevents[j],Ntoys);
      double cor = arraycorD(Nevents[i],Nevents[j],Ntoys);
      cout<<"cor: "<<i<<","<<j<<" = "<<cor<<endl;
     
      hCov->SetBinContent(i+1,j+1,cov);
      hCor->SetBinContent(i+1,j+1,cor);
      //
      if (i!=j){
        hCov->SetBinContent(j+1,i+1,cov);
        hCor->SetBinContent(j+1,i+1,cor);
      }

      // fill diagonal errors
      if (i==j){
        hDiagonalErrors->SetBinContent(i+1,arraymeanD(Nevents[i],Ntoys));
        hDiagonalErrors->SetBinError(i+1,TMath::Sqrt(arrayvarD(Nevents[i],Ntoys)));
      }
    }
  }

  hCor->SetMinimum(-1.11);
  hCor->SetMaximum(1.11);

  drawCor();
  //
  return;
}



///////////////////////////////////////////////////////////////
/*
void SKError::calcAllEff(int ntoy){

  for (int iclass=0; iclass<Nclass; iclass++){

    // 0th index is nominal
    double eff = calcEff(iclass,ntoy,0); 
    cout<<"eff calc: "<<iclass<<","<<ntoy<<": "<<eff<<endl;
    Efficiency[iclass][ntoy] = eff;
  }

  return;
}
*/


///////////////////////////////////////////////////////////////
void SKError::calcAllDelEff(int ntoy){

  for (int iclass=0; iclass<Nclass; iclass++){

    // 0th index is nominal
    double eff = calcDelEff(iclass,ntoy); 
    cout<<"eff calc: "<<iclass<<","<<ntoy<<": "<<eff<<endl;
    DelEfficiency[iclass][ntoy] = eff;
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
//   calcAllEff(ntoy);
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
  Nclass = hEvisNuECCQE->GetNbinsX() + hEvisNuECCOth->GetNbinsX() +  hEvisNuMuCCQE->GetNbinsX() + hEvisNuMuCCOth->GetNbinsX();

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
