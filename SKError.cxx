#ifndef SKERROR_C
#define SKERROR_C

#include "SKError.h"


///////////////////////////////////////////////////////////////
void SKError::initClassTags(){

  // names of classes
  className[0] = "#nu_{e} CC1e [100, 300] MeV";
  className[1] = "#nu_{e} CC1e [300, 700] MeV";
  className[2] = "#nu_{e} CC1e [700, 1250] MeV";
  className[3] = "#nu_{e} CC1e [1250, 2000] MeV";
  className[4] = "#nu_{e} CC1e [2000, 5000] MeV";
  className[5] = "#nu_{e} CC1e [5000, 30000] MeV";
  //
  className[6] = "#nu_{e} CCOth. [0, 1250] MeV";
  className[7] = "#nu_{e} CCOth. [1250, 5000] MeV";
  className[8] = "#nu_{e} CCOth. [5000, 30000] MeV";
  //
  className[9] = "#nu_{#mu} CC1#mu [0, 100] MeV";
  className[10] = "#nu_{#mu} CC1#mu [100, 300] MeV";
  className[11] = "#nu_{#mu} CC1#mu [300, 700] MeV";
  className[12] = "#nu_{#mu} CC1#mu [700, 1250] MeV";
  className[13] = "#nu_{#mu} CC1#mu [1250, 2000] MeV";
  className[14] = "#nu_{#mu} CC1#mu [2000, 5000] MeV";
  className[15] = "#nu_{#mu} CC1#mu [5000, 30000] MeV";
  //
  className[16] = "#nu_{#mu} CCOth. [0, 1250] MeV";
  className[17] = "#nu_{#mu} CCOth. [1250, 5000] MeV";
  className[18] = "#nu_{#mu} CCOth. [5000, 30000] MeV";

  // binnin ames of classes
  classRange[0] = "[100, 300] MeV";
  classRange[1] = "[300, 700] MeV";
  classRange[2] = "[700, 1250] MeV";
  classRange[3] = "[1250, 2000] MeV";
  classRange[4] = "[2000, 5000] MeV";
  classRange[5] = "[5000, 30000] MeV";
  //
  classRange[6] = "[0, 1250] MeV";
  classRange[7] = "[1250, 5000] MeV";
  classRange[8] = "[5000, 30000] MeV";
  //
  classRange[9] = "[0, 100] MeV";
  classRange[10] = "[100, 300] MeV";
  classRange[11] = "[300, 700] MeV";
  classRange[12] = "[700, 1250] MeV";
  classRange[13] = "[1250, 2000] MeV";
  classRange[14] = "[2000, 5000] MeV";
  classRange[15] = "[5000, 30000] MeV";
  //
  classRange[16] = "[0, 1250] MeV";
  classRange[17] = "[1250, 5000] MeV";
  classRange[18] = "[5000, 30000] MeV";

  return;
}

///////////////////////////////////////////////////////////////
void SKError::initTN186Errors(){

  // from TN
  /*
  double shifts[] = {0.555,1.381,-7.092,1.785,4.6263,11.940
                    ,-7.123,5.891,-10.826
                    ,0.87,1.227,0.094,0.136,-0.037,-0.102,-0.157
                    ,6.406, 2.766, -5.67};
  double fits[] = {1.517,2.045,4.380,2.428,7.912,9.272
                  ,10.187,9.096,11.705
                  ,0.376,0.449,0.582,1.194,0.781,1.334,2.485
                  ,5.771,6.109,11.559};
  */

  double nue_ccqe_serr[] = {0.555,1.381,-7.092,1.785,4.6263,11.940};
  double num_ccqe_serr[] = {0.87,1.227,0.094,0.136,-0.037,-0.102,-0.157};
  double nue_ccoth_serr[] = {-7.123,5.891,-10.826};
  double num_ccoth_serr[] = {6.406, 2.766, -5.67};
  //
  double nue_ccqe_ferr[] = {1.517,2.045,4.380,2.428,7.912,9.272};
  double num_ccqe_ferr[] = {0.376,0.449,0.582,1.194,0.781,1.334,2.485};
  double nue_ccoth_ferr[] = {10.187,9.096,11.705};
  double num_ccoth_ferr[] = {5.771,6.109,11.559};
  
  // make histos
  hErrorTN186CCQE[0] = new TH1D("tn186_nueccqe","tn186_nueccqe",7,0,7);
  hErrorTN186CCQE[0]->SetFillColor(kGray);
  hErrorTN186CCQE[0]->SetMarkerStyle(1);
  hErrorTN186CCQE[1] = new TH1D("tn186_numccqe","tn186_numccqe",7,0,7);
  hErrorTN186CCQE[1]->SetFillColor(kGray);
  hErrorTN186CCQE[1]->SetMarkerStyle(1);
  hErrorTN186CCOth[0] = new TH1D("tn186_nueccoth","tn186_nueccoth",3,0,3);
  hErrorTN186CCOth[0]->SetFillColor(kGray);
  hErrorTN186CCOth[0]->SetMarkerStyle(1);
  hErrorTN186CCOth[1] = new TH1D("tn186_numccoth","tn186_numccoth",3,0,3);
  hErrorTN186CCOth[1]->SetFillColor(kGray);
  hErrorTN186CCOth[1]->SetMarkerStyle(1);

  // fill in histos
  // ccqe nue
  for (int i=2; i<=hErrorTN186CCQE[0]->GetNbinsX(); i++){
    hErrorTN186CCQE[0]->SetBinContent(i,nue_ccqe_serr[i-2]);
    hErrorTN186CCQE[0]->SetBinError(i,nue_ccqe_ferr[i-2]);
  }
  // ccqe nue
  for (int i=1; i<=hErrorTN186CCQE[1]->GetNbinsX(); i++){
    hErrorTN186CCQE[1]->SetBinContent(i,num_ccqe_serr[i-1]);
    hErrorTN186CCQE[1]->SetBinError(i,num_ccqe_ferr[i-1]);
  }
  for (int i=1; i<=hErrorTN186CCOth[0]->GetNbinsX(); i++){
    hErrorTN186CCOth[0]->SetBinContent(i,nue_ccoth_serr[i-1]);
    hErrorTN186CCOth[0]->SetBinError(i,nue_ccoth_ferr[i-1]);
  }
  for (int i=1; i<=hErrorTN186CCOth[1]->GetNbinsX(); i++){
    hErrorTN186CCOth[1]->SetBinContent(i,nue_ccoth_serr[i-1]);
    hErrorTN186CCOth[1]->SetBinError(i,nue_ccoth_ferr[i-1]);
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
      vlabel.push_back(new TLatex(offset,xpos,Form("[%1.1f - %1.1f] GeV",evismin,evismax)));
      vlabel.at(i-1)->SetTextAlign(32);
      vlabel.at(i-1)->SetTextSize(0.08);
    }
  }

  return vlabel;
}



///////////////////////////////////////////////////////////////
void SKError::drawErrors1D(bool flgDrawTN186){

  // calc errors first
  calcErrors1D();

  // canvas setup
  TCanvas* cc = new TCanvas("cc","cc",800,600);
  cc->Divide(2,2);
 
  TH1D* herr[] = {hErrorsCCQE[0],hErrorsCCQE[1],hErrorsCCOth[0],hErrorsCCOth[1]};

  // histogram properties
  hErrorsCCQE[0]->SetLineColor(kBlue);;
//  hErrorsCCQE[0]->SetFillColor(kBlue);
  hErrorsCCQE[0]->GetYaxis()->SetTitle("#Delta #varepsilon [%]");
  hErrorsCCQE[0]->SetLineWidth(3);
  hErrorsCCQE[0]->SetMarkerStyle(0);
  hErrorsCCQE[1]->SetLineColor(kRed);
//  hErrorsCCQE[1]->SetFillColor(kRed);
  hErrorsCCQE[1]->GetYaxis()->SetTitle("#Delta #varepsilon [%]");
  hErrorsCCQE[1]->SetLineWidth(3);
  hErrorsCCQE[1]->SetMarkerStyle(0);
  hErrorsCCOth[0]->SetLineColor(kBlue);
//  hErrorsCCOth[0]->SetFillColor(kBlue);
  hErrorsCCOth[0]->GetYaxis()->SetTitle("#Delta #varepsilon [%]");
  hErrorsCCOth[0]->SetLineWidth(3);
  hErrorsCCOth[0]->SetMarkerStyle(0);
  hErrorsCCOth[1]->SetLineColor(kRed);
//  hErrorsCCOth[1]->SetFillColor(kRed);
  hErrorsCCOth[1]->GetYaxis()->SetTitle("#Delta #varepsilon [%]");
  hErrorsCCOth[1]->SetLineWidth(3);
  hErrorsCCOth[1]->SetMarkerStyle(0);

  // set symmetric ranges
  double range_nue_ccqe = TMath::Max( TMath::Abs(hErrorsCCQE[0]->GetMaximum()),
                                      TMath::Abs(hErrorsCCQE[0]->GetMinimum()));
  double range_num_ccqe = TMath::Max( TMath::Abs(hErrorsCCQE[1]->GetMaximum()),
                                      TMath::Abs(hErrorsCCQE[1]->GetMinimum()));
  double range_nue_ccoth = TMath::Max( TMath::Abs(hErrorsCCOth[0]->GetMaximum()),
                                      TMath::Abs(hErrorsCCOth[0]->GetMinimum()));
  double range_num_ccoth = TMath::Max( TMath::Abs(hErrorsCCOth[1]->GetMaximum()),
                                      TMath::Abs(hErrorsCCOth[1]->GetMinimum()));
  double factor = 1.5;
  hErrorsCCQE[0]->SetMinimum(-1*factor*range_nue_ccqe);
  hErrorsCCQE[0]->SetMaximum(1*factor*range_nue_ccqe);
  hErrorsCCQE[1]->SetMinimum(-1*factor*range_num_ccqe);
  hErrorsCCQE[1]->SetMaximum(1*factor*range_num_ccqe);
  hErrorsCCOth[0]->SetMinimum(-1*factor*range_nue_ccoth);
  hErrorsCCOth[0]->SetMaximum(1*factor*range_nue_ccoth);
  hErrorsCCOth[1]->SetMinimum(-1*factor*range_num_ccoth);
  hErrorsCCOth[1]->SetMaximum(1*factor*range_num_ccoth);
  double hmin[] = {-20,-1,-60,-6};
  double hmax[] = {20,1,60,6};
  for (int i=0; i<4; i++){
    herr[i]->SetMinimum(hmin[i]);
    herr[i]->SetMaximum(hmax[i]);
  }
  if (flgDrawTN186){
    double hmin2[] = {-25,-4,-60,-30};
    double hmax2[] = {25,4,60,30};
    for (int i=0; i<4; i++){
      herr[i]->SetMinimum(hmin2[i]);
      herr[i]->SetMaximum(hmax2[i]);
    }
  }

  // get bin labels
  vector<TLatex*> vlabel_ccqe = getBinLabels(hEvisNuMuCCQE);
  vector<TLatex*> vlabel_ccoth = getBinLabels(hEvisNuECCOth);
  // set bin labels
  for (int ibin=1; ibin<=hErrorsCCQE[0]->GetNbinsX(); ibin++){
    hErrorsCCQE[0]->GetXaxis()->SetBinLabel(ibin,vlabel_ccqe.at(ibin-1)->GetTitle());
    hErrorsCCQE[1]->GetXaxis()->SetBinLabel(ibin,vlabel_ccqe.at(ibin-1)->GetTitle());
  }
  for (int ibin=1; ibin<=hErrorsCCOth[0]->GetNbinsX(); ibin++){
    hErrorsCCOth[0]->GetXaxis()->SetBinLabel(ibin,vlabel_ccoth.at(ibin-1)->GetTitle());
    hErrorsCCOth[1]->GetXaxis()->SetBinLabel(ibin,vlabel_ccoth.at(ibin-1)->GetTitle());
  }
  float labelsize = 0.06;
  hErrorsCCQE[0]->GetXaxis()->SetLabelSize(labelsize);
  hErrorsCCQE[1]->GetXaxis()->SetLabelSize(labelsize);
  labelsize = 0.08;
  hErrorsCCOth[0]->GetXaxis()->SetLabelSize(labelsize);
  hErrorsCCOth[1]->GetXaxis()->SetLabelSize(labelsize);

  // draw
  TF1* zeroValue = new TF1("zero","0",0,10000);
  zeroValue->SetLineStyle(2);
  zeroValue->SetLineWidth(2);
  cc->cd(1);
  hErrorsCCQE[0]->Draw("EH");
  zeroValue->Draw("same");
  cc->cd(2);
  hErrorsCCQE[1]->Draw("EH");
  zeroValue->Draw("same");
  cc->cd(3);
  hErrorsCCOth[0]->Draw("EH");
  zeroValue->Draw("same");
  cc->cd(4);
  hErrorsCCOth[1]->Draw("EH");
  zeroValue->Draw("same");
  if (flgDrawTN186){
    cc->cd(1);
    hErrorTN186CCQE[0]->Draw("E2same");
    hErrorsCCQE[0]->Draw("sameEH");
    zeroValue->Draw("same");
    cc->cd(2);
    hErrorTN186CCQE[1]->Draw("E2same");
    hErrorsCCQE[1]->Draw("sameEH");
    zeroValue->Draw("same");
    cc->cd(3);
    hErrorTN186CCOth[0]->Draw("E2same");
    hErrorsCCOth[0]->Draw("sameEH");
    zeroValue->Draw("same");
    cc->cd(4);
    hErrorTN186CCOth[1]->Draw("E2same");
    hErrorsCCOth[1]->Draw("sameEH");
    zeroValue->Draw("same");
  }
 

  //
  return;

}



///////////////////////////////////////////////////////////////
void SKError::calcErrors1D(){

  // identify bins
  int bins_nue_ccqe[] = {1,2,3,4,5,6};
  int bins_nue_ccoth[] = {7,8,9};
  int bins_num_ccqe[] = {10,11,12,13,14,15,16};
  int bins_num_ccoth[] = {17,18,19};

  // convert to percent
  float factor = 100.0;

  // ccqe nue
  for (int i=hErrorsCCQE[0]->GetNbinsX(); i>1; i--){
    int nbin = bins_nue_ccqe[i-2];
    hErrorsCCQE[0]->SetBinContent(i,factor*hErrors->GetBinContent(nbin));
    hErrorsCCQE[0]->SetBinError(i,factor*hErrors->GetBinError(nbin));
  }

  // ccqe numu
  for (int i=hErrorsCCQE[1]->GetNbinsX(); i>0; i--){
    // fill numu ccqe bins
    int nbin = bins_num_ccqe[i-1];
    hErrorsCCQE[1]->SetBinContent(i,factor*hErrors->GetBinContent(nbin));
    hErrorsCCQE[1]->SetBinError(i,factor*hErrors->GetBinError(nbin));
  }

  // ccother nue
  for (int i=hErrorsCCOth[0]->GetNbinsX(); i>0; i--){
    // fill numu ccqe bins
    int nbin = bins_nue_ccoth[i-1];
    hErrorsCCOth[0]->SetBinContent(i,factor*hErrors->GetBinContent(nbin));
    hErrorsCCOth[0]->SetBinError(i,factor*hErrors->GetBinError(nbin));
  }

  // ccother numu
  for (int i=hErrorsCCOth[1]->GetNbinsX(); i>0; i--){
    // fill numu ccqe bins
    int nbin = bins_num_ccoth[i-1];
    hErrorsCCOth[1]->SetBinContent(i,factor*hErrors->GetBinContent(nbin));
    hErrorsCCOth[1]->SetBinError(i,factor*hErrors->GetBinError(nbin));
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
void SKError::marginalize(int ntoy){

  // precondition: Nevents arrays filled over marginal points
  for (int iclass=0; iclass<Nclass; iclass++){
    DelEfficiency[iclass][ntoy] = calcMargEff(iclass);
  }

  //
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
void SKError::printAllEffDist(const char* plotdir, int nutype, const char* tag){
     TString stag = tag;

     // nutype = 1 (electron) or 2 (muon)

     // set up list of classes
     int iclass[10];
     int nclasses = 9;
     if (nutype == 1){
        int start = 0;
        for ( int i=0; i<nclasses; i++ ) {
          iclass[i] = start + i;
        }
     }
     else {
        nclasses = 10;
        int start = 9;
        for ( int i=0; i<nclasses; i++ ) {
          iclass[i] = start + i;
        }
     }

     // set up canvas
     TCanvas* cc = new TCanvas("cc","cc",900,600);
     cc->Divide(4,3);

     // set up pads 
     int ipad[10];
     if (nutype == 1){
       ipad[0] = 1; 
       ipad[1] = 2; 
       ipad[2] = 3; 
       ipad[3] = 5; 
       ipad[4] = 6; 
       ipad[5] = 7; 
       ipad[6] = 9; 
       ipad[7] = 10; 
       ipad[8] = 11; 
     }
     else {
       ipad[0] = 1; 
       ipad[1] = 2; 
       ipad[2] = 3; 
       ipad[3] = 4; 
       ipad[4] = 5; 
       ipad[5] = 6; 
       ipad[6] = 7; 
       ipad[7] = 9; 
       ipad[8] = 10; 
       ipad[9] = 11; 
     }


     // directory to dump plots
     TString pdir = plotdir;

     // print them
     for (int i=0; i<nclasses; i++){
       // got to right pad
       cc->cd(ipad[i]);
       // draw distribution
       drawEffDist(iclass[i]);
     }


     TString pname = pdir.Data();
     pname.Append(Form("eff_dist_nutype%d",nutype));
     pname.Append(stag.Data());
     pname.Append(".png");
     cc->Print(pname.Data());

     return;
}



///////////////////////////////////////////////////////////////
void SKError::printEachEffDist(const char* plotdir){
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
  TString sector[] = {"CC1e","CCOth.","CC1#mu","CCOth."};
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
  double eff = Nevents[nclass][ntoy] / NeventsNominal[nclass];
  cout<<"eff: "<<eff<<endl;
  double eff0 = NeventsNominalCore[nclass] / NeventsNominal[nclass];
  cout<<"eff0: "<<eff0<<endl;

  // get this efficiency  
//  double eff = Nevents[nclass][ntoy] / NeventsTotal[nclass][0];
//  cout<<"eff: "<<eff<<endl;
//  double eff0 = Nevents[nclass][0] / NeventsTotal[nclass][0];
  cout<<"eff0: "<<eff0<<endl;
  
  if (eff0 == 0) return 0.;

  return ((eff/eff0) - 1.); 

}




///////////////////////////////////////////////////////////////
// calcutlate marginalized efficiency
// here the modified eff. is calc w.r.t. modified total
double SKError::calcMargEff(int nclass){

  // fill array of efficiencies
  double effs[NTOYS];
  for (int i=0; i<Nmarginal; i++){
    double eff = Nevents[nclass][i] / NeventsTotal[nclass][i];
    double eff0 = NeventsNominalCore[nclass] / NeventsNominal[nclass];
    effs[i] = ((eff/eff0) - 1.); 
  }

  // return average
  return arraymeanD(effs,Nmarginal);
};



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
  double eff0 = NeventsNominalCore[nclass] / NeventsNominal[nclass];
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
      hErrors->SetBinContent(i+1,arraymeanD(DelEfficiency[i],Ntoys));
      hErrors->SetBinError(i+1,TMath::Sqrt(arrayvarD(DelEfficiency[i],Ntoys)));
      double cov = arraycovD(DelEfficiency[i],DelEfficiency[i],Ntoys);
      double cor = arraycorD(DelEfficiency[i],DelEfficiency[i],Ntoys);
      hCov->SetBinContent(i+1,i+1,cov);
      hCor->SetBinContent(i+1,i+1,cor);
    }

    // NSK errors
    else if ( vartype==1 ){
      hErrors->SetBinContent(i+1,arraymeanD(Nevents[i],Ntoys));
      hErrors->SetBinError(i+1,TMath::Sqrt(arrayvarD(Nevents[i],Ntoys)));       
      double cov = arraycovD(Nevents[i],Nevents[i],Ntoys);
      double cor = arraycorD(Nevents[i],Nevents[i],Ntoys);
      hCov->SetBinContent(i+1,i+1,cov);
      hCor->SetBinContent(i+1,i+1,cor);
    }

    // fill off-diagonal errors
    for (int j=0; j<Nclass; j++){
     
      // symmetry
      if (j<i) continue;
         
      // off-diagonal elements
      if ( vartype==0 ) {
        double cov = arraycovD(DelEfficiency[i],DelEfficiency[j],Ntoys);
        double cor = arraycorD(DelEfficiency[i],DelEfficiency[j],Ntoys);
        if (j==i){
          cor = 1.;
        }
        hCov->SetBinContent(i+1,j+1,cov);
        hCor->SetBinContent(i+1,j+1,cor);
        hCov->SetBinContent(j+1,i+1,cov);
        hCor->SetBinContent(j+1,i+1,cor);
      }
      else if ( vartype==1 ){
        double cov = arraycovD(Nevents[i],Nevents[j],Ntoys);
        double cor = arraycorD(Nevents[i],Nevents[j],Ntoys);
        if (j==i){
          cor = 1.;
        }
        hCov->SetBinContent(i+1,j+1,cov);
        hCor->SetBinContent(i+1,j+1,cor);
        hCov->SetBinContent(j+1,i+1,cov);
        hCor->SetBinContent(j+1,i+1,cor);
      }
    }
  }

  //scales
  hCor->SetMinimum(-1.05);
  hCor->SetMaximum(1.00);

  // draw correlations
  drawCor();

  return;
}



///////////////////////////////////////////////////////////////
void SKError::calcAllEff(int ntoy){

  // find efficiency for each class
  for (int iclass=0; iclass<Nclass; iclass++){

    double eff = 0.;
    if (effDefinition==0){
      eff = calcDelEff(iclass,ntoy); 
      DelEfficiency[iclass][ntoy] = eff;
    }

    // use nominal value (instead of alpha-modified, see TN-186)
    // in the denominator:
//    else if (effDefinition==1){
//      eff = calcEff(iclass,ntoy); 
//      DelEfficiency[iclass][ntoy] = eff;
//    }

  }

  return;
}



///////////////////////////////////////////////////////////////
void SKError::addNominal(){

   int index = 0;
   for (int ih=0; ih<NEvisHistos; ih++){
     for (int ibin=1; ibin<=hEvisCore[ih]->GetNbinsX(); ibin++ ){
       NeventsNominalCore[index] = hEvisCore[ih]->GetBinContent(ibin);
       NeventsNominal[index] = hEvisBase[ih]->GetBinContent(ibin);
       index++;
     }
   }
 
   /*  
   // save histos in arrays
   int index = 0;
   //
   for (int ibin=1; ibin<=hEvisNuECCQE->GetNbinsX(); ibin++){
     NeventsNominalCore[index] = hEvisNuECCQE->GetBinContent(ibin);
     index++;
   }
   //
   for (int ibin=1; ibin<=hEvisNuECCOth->GetNbinsX(); ibin++){
     NeventsNominalCore[index] = hEvisNuECCOth->GetBinContent(ibin);
     index++;
   }
   //
   for (int ibin=1; ibin<=hEvisNuMuCCQE->GetNbinsX(); ibin++){
     NeventsNominalCore[index] = hEvisNuMuCCQE->GetBinContent(ibin);
     index++;
   }
   //;
   for (int ibin=1; ibin<=hEvisNuMuCCOth->GetNbinsX(); ibin++){
     NeventsNominalCore[index] = hEvisNuMuCCOth->GetBinContent(ibin);
     index++;
   }

   // save histos in TOTAL arrays
   index = 0;
   //
   for (int ibin=1; ibin<=hEvisNuECCQE->GetNbinsX(); ibin++){
     NeventsNominal[index] = hEvisNuECCQETot->GetBinContent(ibin);
     index++;
   }
   //
   for (int ibin=1; ibin<=hEvisNuECCOth->GetNbinsX(); ibin++){
     NeventsNominal[index] = hEvisNuECCOthTot->GetBinContent(ibin);
     index++;
   }
   //
   for (int ibin=1; ibin<=hEvisNuMuCCQE->GetNbinsX(); ibin++){
     NeventsNominal[index] = hEvisNuMuCCQETot->GetBinContent(ibin);
     index++;
   }
   //
   for (int ibin=1; ibin<=hEvisNuMuCCOth->GetNbinsX(); ibin++){
     NeventsNominal[index] = hEvisNuMuCCOthTot->GetBinContent(ibin);
     index++;
   }
   
   // calculate all of the epsilon values!
//   calcAllEff(ntoy);
  */

 // 
 return;


}



///////////////////////////////////////////////////////////////
void SKError::fillNeventArrays(int ntoy){

   int index = 0;
   for (int ih=0; ih<NEvisHistos; ih++){
     for (int ibin=1; ibin<=hEvisCore[ih]->GetNbinsX(); ibin++ ){
       cout<<"total events index: "<<index<<"content "<<hEvisBase[ih]->GetBinContent(ibin)<<endl;
       Nevents[index][ntoy] = hEvisCore[ih]->GetBinContent(ibin);
       NeventsTotal[index][ntoy] = hEvisBase[ih]->GetBinContent(ibin);
       index++;
     }
   }

   /*  
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
   */
   return;
}



///////////////////////////////////////////////////////////////
void SKError::addToy(int ntoy){

   // fill all event arrays
   fillNeventArrays(ntoy);

   // calculate all of the epsilon values!
   calcAllEff(ntoy);

 // 
 return;

}


///////////////////////////////////////////////////////////////
void SKError::resetHistos(){

  for (int ih=0; ih<NEvisHistos; ih++){
    cout<<"SKError::restHistos: resetting histograms: "<<ih<<endl;
    hEvisCore[ih]->Reset();
    hEvisBase[ih]->Reset();
  }

  //
//  hEvisNuECCQE->Reset();
//  hEvisNuMuCCQE->Reset();
//  hEvisNuECCOth->Reset();
//  hEvisNuMuCCOth->Reset();
  //
//  hEvisNuECCQETot->Reset();
//  hEvisNuMuCCQETot->Reset();
//  hEvisNuECCOthTot->Reset();
//  hEvisNuMuCCOthTot->Reset();



  return;
}



///////////////////////////////////////////////////////////////
int SKError::addEventCore(int nclass, double evis, double weight, int coretype){
 
    //
    if (nclass==1){
      return hEvisNuECCQE->Fill(evis,weight);
    }
    if (nclass==2){
      return hEvisNuMuCCQE->Fill(evis,weight);
    }
    if (nclass==3){
      return hEvisNuECCOth->Fill(evis,weight);
    }
    if (nclass==4){
      return hEvisNuMuCCOth->Fill(evis,weight);
    }
    if (flgNC){
      if (nclass==5 && coretype==1){
        hEvisNuENCPi0->Fill(evis,weight);
      }
      else if (nclass==5 && coretype==2){
        hEvisNuMuNCPi0->Fill(evis,weight);
      }
      else if (nclass==6 && coretype==1){
        hEvisNuENCPi->Fill(evis,weight);
      }
      else if (nclass==6 && coretype==2){
        hEvisNuMuNCPi->Fill(evis,weight);
      }
    }

  return 0;
}



///////////////////////////////////////////////////////////////
int SKError::addEventBase(int nclass, double evis, double weight){

  
    //
    if (nclass==1){
      return hEvisNuECCQETot->Fill(evis,weight);
    }
    if (nclass==2){
      return hEvisNuMuCCQETot->Fill(evis,weight);
    }
    if (nclass==3){
      return hEvisNuECCOthTot->Fill(evis,weight);
    }
    if (nclass==4){
      return hEvisNuMuCCOthTot->Fill(evis,weight);
    }
    if (flgNC){
      if (nclass==5){
        hEvisNuENCPi0Tot->Fill(evis,weight);
        hEvisNuMuNCPi0Tot->Fill(evis,weight);
      }
      if (nclass==6){
        hEvisNuENCPiTot->Fill(evis,weight);
        hEvisNuMuNCPiTot->Fill(evis,weight);
      }
    }

  return 0;
}

///////////////////////////////////////////////////////////////
int SKError::addEvent(int nclass, double evis, double weight, bool flgtot){

  
  // passes core
  if (!flgtot){
    //
    if (nclass==1){
      return hEvisNuECCQE->Fill(evis,weight);
    }
    if (nclass==2){
      return hEvisNuMuCCQE->Fill(evis,weight);
    }
    if (nclass==3){
      return hEvisNuECCOth->Fill(evis,weight);
    }
    if (nclass==4){
      return hEvisNuMuCCOth->Fill(evis,weight);
    }
  }

  // total
  else{
    //
    if (nclass==1){
      return hEvisNuECCQETot->Fill(evis,weight);
    }
    if (nclass==2){
      return hEvisNuMuCCQETot->Fill(evis,weight);
    }
    if (nclass==3){
      return hEvisNuECCOthTot->Fill(evis,weight);
    }
    if (nclass==4){
      return hEvisNuMuCCOthTot->Fill(evis,weight);
    }
  }

  return 0;
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
  //  6 -> NC pi

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

  // NC pi0
  if (component==4 && TMath::Abs(mode)>=30 && evis>100. && nsubev<=2
       && towall>towall_nuecut && wall>wall_nuecut ){
    return 5;
  }

  // NC pi
  if (component==5 && TMath::Abs(mode)>=30 && evis>30. && nsubev<=2
       && towall>towall_numucut && wall>wall_numucut ){
    return 6;
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
   double xmax = arraymaxD(DelEfficiency[iclass],Ntoys);
   double width = (xmax-xmin);
   xmin = xmin-(width/2.);
   xmax = xmax+(width/2.);

   // make histogram
   int nbinsx = 30;
   hdist=new TH1D("hdist","hdist",nbinsx,xmin,xmax);
   for (int i=0; i<Ntoys; i++){
     hdist->Fill(DelEfficiency[iclass][i]);
   }

   // setup histogram
   hdist->SetLineWidth(3);
   hdist->SetTitle(className[iclass].Data());
   hdist->SetTitleSize(0.25);
   hdist->SetLineColor(kBlue);
   hdist->GetXaxis()->SetTitle("#Delta #varepsilon");
//   hdist->GetYaxis()->SetTitle("# of throws");
   hdist->GetYaxis()->SetTitleSize(0.065);
   hdist->GetXaxis()->SetTitleSize(0.065);
   hdist->GetXaxis()->SetTitleOffset(0.7);
   hdist->GetYaxis()->SetTitleOffset(0.7);
   hdist->GetXaxis()->SetNdivisions(5);
   hdist->GetXaxis()->SetLabelSize(0.06);
   hdist->GetYaxis()->SetLabelSize(0.06);
   hdist->GetYaxis()->SetNdivisions(5);
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

  // histo binning )
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

  hEvisCore[0]=hEvisNuECCQE;
  hEvisCore[1]=hEvisNuECCOth;
  hEvisCore[2]=hEvisNuMuCCQE;
  hEvisCore[3]=hEvisNuMuCCOth;
  //
  hEvisBase[0]=hEvisNuECCQETot;
  hEvisBase[1]=hEvisNuECCOthTot;
  hEvisBase[2]=hEvisNuMuCCQETot;
  hEvisBase[3]=hEvisNuMuCCOthTot;

  NEvisHistos = 4;
  }

  if (ibinning==0){
  // histo binning )
  int    NbinsNuECCQE = 2;
  double BinsNuECCQE[] = {100,1250.,30000.};
  int    NbinsNuECCOth = 2;
  double BinsNuECCOth[] = {100.,1250.,30000.};
  int    NbinsNuMuCCQE = 3;
  double BinsNuMuCCQE[] = {0.,1250.,5000.,30000.};
  int    NbinsNuMuCCOth = 2;
  double BinsNuMuCCOth[] = {0,1250.,30000.};
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

  hEvisCore[0]=hEvisNuECCQE;
  hEvisCore[1]=hEvisNuECCOth;
  hEvisCore[2]=hEvisNuMuCCQE;
  hEvisCore[3]=hEvisNuMuCCOth;
  //
  hEvisBase[0]=hEvisNuECCQETot;
  hEvisBase[1]=hEvisNuECCOthTot;
  hEvisBase[2]=hEvisNuMuCCQETot;
  hEvisBase[3]=hEvisNuMuCCOthTot;

  NEvisHistos = 4;
  }



  if (ibinning==2){
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
  hEvisNuENCPi0 = new TH1D("hnuencpi0","hnuencpi0",1,100,1250.); 
  hEvisNuENCPi = new TH1D("hnuencpi","hnuencpi",1,100,1250.); 
  hEvisNuMuNCPi0 = new TH1D("hnumncpi0","hnumncpi0",1,30,3000.); 
  hEvisNuMuNCPi = new TH1D("hnumncpi","hnumncpi",1,30,3000.); 

  // totals
  hEvisNuECCQETot = new TH1D("hnueccqetot","hnueccqetot",NbinsNuECCQE,BinsNuECCQE); 
  hEvisNuECCOthTot = new TH1D("hnueccothtot","hnueccothtot",NbinsNuECCOth,BinsNuECCOth); 
  hEvisNuMuCCQETot = new TH1D("hnumuccqetot","hnumuccqetot",NbinsNuMuCCQE,BinsNuMuCCQE); 
  hEvisNuMuCCOthTot = new TH1D("hnumuccothtot","hnumuccothtot",NbinsNuMuCCOth,BinsNuMuCCOth); 
  hEvisNuENCPi0Tot = new TH1D("hnuencpi0tot","hnuencpi0tot",1,30,30000.); 
  hEvisNuENCPiTot = new TH1D("hnuencpitot","hnuencptoti",1,30,30000.); 
  hEvisNuMuNCPi0Tot = new TH1D("hnumncpi0tot","hnumncpi0tot",1,30,30000.); 
  hEvisNuMuNCPiTot = new TH1D("hnumncpitot","hnumncpitot",1,30,30000.); 

  hEvisCore[0]=hEvisNuECCQE;
  hEvisCore[1]=hEvisNuECCOth;
  hEvisCore[2]=hEvisNuMuCCQE;
  hEvisCore[3]=hEvisNuMuCCOth;
  hEvisCore[4]=hEvisNuENCPi0;
  hEvisCore[5]=hEvisNuENCPi;
  hEvisCore[6]=hEvisNuMuNCPi0;
  hEvisCore[7]=hEvisNuMuNCPi;
  //
  hEvisBase[0]=hEvisNuECCQETot;
  hEvisBase[1]=hEvisNuECCOthTot;
  hEvisBase[2]=hEvisNuMuCCQETot;
  hEvisBase[3]=hEvisNuMuCCOthTot;
  hEvisBase[4]=hEvisNuENCPi0Tot;
  hEvisBase[5]=hEvisNuENCPiTot;
  hEvisBase[6]=hEvisNuMuNCPi0Tot;
  hEvisBase[7]=hEvisNuMuNCPiTot;

  flgNC = true;
  NEvisHistos = 8;
  }

  // count all bins
  Nclass=0;
  for (int ih=0; ih<NEvisHistos; ih++){
   Nclass += hEvisCore[ih]->GetNbinsX(); 
  }
  // count all bins
//  Nclass = hEvisNuECCQE->GetNbinsX() + hEvisNuECCOth->GetNbinsX()
//           +  hEvisNuMuCCQE->GetNbinsX() + hEvisNuMuCCOth->GetNbinsX();

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
  hErrors = new TH1D("hdiag","hdiag",Nclass,0,Nclass);
  int nbins_ccqe = TMath::Max(hEvisNuECCQE->GetNbinsX(),hEvisNuMuCCQE->GetNbinsX());
  int nbins_ccoth = TMath::Max(hEvisNuECCOth->GetNbinsX(),hEvisNuMuCCOth->GetNbinsX());
  for ( int i=0; i<2; i++ ) {
    hErrorsCCQE[i] = new TH1D("hdiag_ccqe","hdiag_ccqe",nbins_ccqe,0,nbins_ccqe);
    hErrorsCCOth[i] = new TH1D("hdiag_ccoth","hdiag_ccoth",nbins_ccoth,0,nbins_ccoth);
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
  initTN186Errors();
  initClassTags();

  for (int iclass=0; iclass<NCLASSES; iclass++){
    NeventsNominal[iclass]=0.;
    NeventsNominalCore[iclass]=0.;
    DelEffShiftError[iclass]=0.;
    absDelEffShiftError[iclass]=0.;
    DelEffFitError[iclass]=0.;
    for (int jtoy=0; jtoy<NTOYS; jtoy++){
      Nevents[iclass][jtoy]=0.;
      NeventsTotal[iclass][jtoy]=0.;
      DelEfficiency[iclass][jtoy]=0.;
    }
  }

  //
  evisClass[1][1]=0;
  evisClass[1][2]=1;
  evisClass[1][3]=2;
  evisClass[1][4]=3;
  evisClass[1][5]=4;
  evisClass[1][6]=5;
  //
  evisClass[3][1]=6;
  evisClass[3][2]=7;
  evisClass[3][3]=8;
  //
  evisClass[2][1]=9;
  evisClass[2][2]=10;
  evisClass[2][3]=11;
  evisClass[2][4]=12;
  evisClass[2][5]=13;
  evisClass[2][6]=14;
  evisClass[2][7]=15;
  //
  evisClass[4][1]=16;
  evisClass[4][2]=17;
  evisClass[4][3]=18;
  //


  // 0 -> 100% defined using flux, xsec, and norm parameters
  // 1 -> 100% defined using nominal value
  effDefinition = 0;

  flgNC = false;
}




///////////////////////////////////////////////////////////////
TH1D* SKError::makeHisto(const char* varname, const char* hname){

  TString var = varname;
  int nbins = 150;
  double xmin = -1000;
  double xmax = 1000;
  TString xtitle = "";

  if (!var.CompareTo("fqpidpar")){
    nbins = 150;
    xmin = -4500;
    xmax = 4500;
    xtitle.Append("e/#mu PID");
  }
  if (!var.CompareTo("fqpi0par")){
    nbins = 100;
    xmin = -1500;
    xmax = 1500;
    xtitle.Append("e/#pi_{0} PID");
  }
  if (!var.CompareTo("fqpippar")){
    nbins = 150;
    xmin = -1500;
    xmax = 1500;
    xtitle.Append("#mu/#pi PID");
  }
  if (!var.CompareTo("fqrcpar")){
    nbins = 150;
    xmin = -400;
    xmax = 1000;
    xtitle.Append("RC parameter");
  }
  if (!var.CompareTo("pmomv")){
    nbins = 100;
    xmin = 0;
    xmax = 30000;
    xtitle.Append("E_{#nu}");
  }
  if (!var.CompareTo("nse")){
    nbins = 10;
    xmin = 0;
    xmax = 10;
    xtitle.Append("N_{subev}");
  }
  if (!var.CompareTo("nring")){
    nbins = 10;
    xmin = 0;
    xmax = 10;
    xtitle.Append("N_{rings}");
  }
  if (!var.CompareTo("mode")){
    nbins = 80;
    xmin = 0;
    xmax = 40;
    xtitle.Append("NEUT mode");
  }
  if (!var.CompareTo("wall")){
    nbins = 100;
    xmin = 0;
    xmax = 1500;
    xtitle.Append("wall");
  }

  TH1D* h = new TH1D(hname,hname,nbins,xmin,xmax);
  h->GetXaxis()->SetTitle(xtitle.Data());
  return h;
}



///////////////////////////////////////////////////////////////
int SKError::isCore(int nclass, float fqpidpar, float fqpi0par, float fqpippar,float fqrcpar){

 // see if it passes core cuts
  // classes: 1 -> nu e CCQE
  //          2 -> nu mu CCQE
  //          3 -> nu e CCOth
  //          4 -> nu mu CCOth
  //          5 -> NC pi0
  if (nclass==1){ //< require single ring electron and CCQE
    if (fqpidpar>=0. && fqpi0par<=0. && fqrcpar<=0.){//< e-like, not pi0, 1R-like
      return 1;
    }
    else{
      return 0;
    }
  }


  // muon CCQE
  if (nclass==2){
    if (fqpidpar<=0 && fqpippar<=0. && fqrcpar<=0.){
      return 1; // mu-lik and not pip and 1R-like so is core
    }
    else{
      return 0;
    }
  }


  // electron CCOth
  if (nclass==3){
    if (fqpidpar>=0 && fqpi0par<=0. && fqrcpar<=0.){
      return 1; // e-like and not pi0 and 1R-like
    }
    else{
      return 0;
    }
  }


  // muon CCOth
  if (nclass==4){
    if (fqpidpar<=0 && fqpippar<=0. && fqrcpar<=0.){
      return 1; //< mu-lik and not pip and 1R-like
    }
    else{
      return 0.;
    }
  }


  // NC
  if (nclass==5){
    if (fqpidpar>=0 && fqpi0par<=0. && fqrcpar<=0.){
      return 1; //< e-like and not pi0-like
    }
    else{
      return 0.;
    }
  }

  //
  return 0;


}



///////////////////////////////////////////////////////////////
int SKError::getClassEvis(int nclass, int nevisbin){
  return evisClass[nclass][nevisbin];
}


///////////////////////////////////////////////////////////////
void SKError::drawVariable(mcLargeArray* events, const char* varname, int nclass, int flgcore){

  // get name of var
  TString var = varname;

  // setup histogram binning
  if (hVariable1!=NULL){
    hVariable1->Delete();
  }
  if (hVariable2!=NULL){
    hVariable2->Delete();
  }
  hVariable2 = makeHisto(varname,"vhvar2");
  hVariable1 = makeHisto(varname,"vhvar1");

  // loop over MC events
  for (int iev=0; iev<events->nsize; iev++){

    // get the class
    int iclass = getClassMC(events->vnutype[iev], events->vmode[iev], events->vcomponent[iev],
                 events->vfqemom[iev], events->vfqnsubev[iev], events->vfqtowall[iev], events->vfqwall[iev]);

    if (iclass==0) continue;
    // get the erec bin
    int ebin = addEvent(iclass,events->vfqemom[iev], 1.,false);
//    int ebin = addEventCore(iclass,events->vfqemom[iev], 1.);
    
    int evisclass = getClassEvis(iclass,ebin);
   

    if (evisclass!=nclass){ 
      continue;
    }

    if (flgcore==2){
      double ww = events->vweight[iev];
      if  (!var.CompareTo("fqpidpar")){
        hVariable2->Fill(events->vattribute[iev][0],ww);
      }
      if  (!var.CompareTo("fqpi0par")){
        hVariable2->Fill(events->vattribute[iev][1],ww);
      }
      if  (!var.CompareTo("fqpippar")){
        hVariable2->Fill(events->vattribute[iev][2],ww);
      }
      if  (!var.CompareTo("fqrcpar")){
        hVariable2->Fill(events->vattribute[iev][3],ww);
      }
      if  (!var.CompareTo("pmomv")){
        hVariable2->Fill(events->vpmomv[iev],ww);
      }
      if  (!var.CompareTo("nse")){
        hVariable2->Fill(events->vfqnsubev[iev],ww);
      }
      if  (!var.CompareTo("mode")){
        hVariable2->Fill(events->vmode[iev],ww);
      }
      if  (!var.CompareTo("wall")){
        hVariable2->Fill(events->vfqwall[iev],ww);
      }
      if  (!var.CompareTo("nring")){
        hVariable2->Fill(events->vfqnring[iev],ww);  
      }
    }
    

    if (flgcore){
      if (!isCore(iclass,events->vattribute[iev][0]
                        ,events->vattribute[iev][1]
                        ,events->vattribute[iev][2]
                        ,events->vattribute[iev][3])){
        continue;
      }
    }

    double ww = events->vweight[iev];
    if  (!var.CompareTo("fqpidpar")){
      hVariable1->Fill(events->vattribute[iev][0],ww);
    }
    if  (!var.CompareTo("fqpi0par")){
      hVariable1->Fill(events->vattribute[iev][1],ww);
    }
    if  (!var.CompareTo("fqpippar")){
      hVariable1->Fill(events->vattribute[iev][2],ww);
    }
    if  (!var.CompareTo("fqrcpar")){
      hVariable1->Fill(events->vattribute[iev][3],ww);
    }
    if  (!var.CompareTo("pmomv")){
      hVariable1->Fill(events->vpmomv[iev],ww);
    }
    if  (!var.CompareTo("nse")){
      hVariable1->Fill(events->vfqnsubev[iev],ww);
    }
    if  (!var.CompareTo("mode")){
      hVariable1->Fill(events->vmode[iev],ww);
    }
    if  (!var.CompareTo("wall")){
      hVariable1->Fill(events->vfqwall[iev],ww);
    }
    if  (!var.CompareTo("nring")){
      hVariable1->Fill(events->vfqnring[iev],ww);
    }
  
  }

  hVariable1->Draw("");
  if (flgcore==2) {
    hVariable2->SetLineColor(kRed);
    hVariable2->Draw("h");
    hVariable1->Draw("sameh");

  }


};



#endif
