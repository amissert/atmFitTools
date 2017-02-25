#ifndef HISOTARRAYFV_H
#define HISOTARRAYFV_H

#include "modHistoArrayFV.h"



using namespace std;

////////////////////////////////////////////////////////////////////////
// Save the overall uncertainties to be used in moreUncertainties.cxx
void modHistoArrayFV::saveSummary(const char* fname){

  // set up file
  TString outfilename = fname;
  TFile* outfile = new TFile(outfilename.Data(),"RECREATE");
  outfile->cd();
  cout<<"Writing summary file: "<<outfile->GetName()<<endl;

  // write histograms
  
  // uncertainty maps
  FVUncMap->Write();
  FVShiftMap->Write();
  FVFitMap->Write();
  //
  FVUncMapCCQE->Write();
  FVUncMapCCnQE->Write();
  FVUncMapCCWrong->Write();
  FVUncMapNC->Write();
  //
  FVShiftMapCCQE->Write();
  FVShiftMapCCnQE->Write();
  FVShiftMapCCWrong->Write();
  FVShiftMapNC->Write();
  //
  FVFitMapCCQE->Write();
  FVFitMapCCnQE->Write();
  FVFitMapCCWrong->Write();
  FVFitMapNC->Write();
  // nominal values
  hFV[0]->Write();
  hFVCCQE[0]->Write();
  hFVCCnQE[0]->Write();
  hFVCCWrong[0]->Write();
  hFVNC[0]->Write();
  // throw histograms
  for (int fvbin=0; fvbin<hFV[0]->GetNumberOfBins(); fvbin++){
   hNevents[fvbin]->Write();
   hNeventsCCQE[fvbin]->Write();
   hNeventsCCnQE[fvbin]->Write();
   hNeventsCCWrong[fvbin]->Write();
   hNeventsNC[fvbin]->Write();
  }

  /*
  for (int ibin=0; ibin<hFV[0]->GetNumberOfBins(); ibin++){
    for (int ebin=0; ebin<binUnc[ibin]->GetNbinsX(); ebin++){
      float binc =  binUnc[ibin]->GetBinContent(ebin);
      if (binc>0.){
        // get "shift" error
        float shift = TMath::Abs(binc - histos[0][ibin]->GetBinContent(ebin));
        // total fractional error is shift + uncertainty /nev
        binUnc[ibin]->SetBinContent(ebin,(binUnc[ibin]->GetBinError(ebin)+shift)/binc);
        binUnc[ibin]->SetBinError(ebin,0.);
      }
    }
    binUnc[ibin]->Write();
//    histos[0][ibin]->Write();
  }
  */
 
  cout<<"Closing summary file: "<<outfile->GetName()<<endl;
  outfile->Close();
  cout<<"Closed!: "<<outfile->GetName()<<endl;
  
  //
  return;
}



///////////////////////////////////////
// save histograms and quit
void modHistoArrayFV::saveClose(){

//  cout<<"Opening array file: "<<fout->GetName()<<endl;
//  fout->Write();
//  cout<<"Closing array file: "<<fout->GetName()<<endl;
  fout->Close();

  return;
}



///////////////////////////////////////////////////////////
// Calculate summary statistics
// This includes the covariance and correlations between the bins
// as well as the overall uncertainty, which is taken as the RMS
// of the number of events in each bin
void modHistoArrayFV::calcSummary(){

  /////////////////////////////////////////////////////
  // find mean for each erec bin
  const int nfvbins = hFV[0]->GetNumberOfBins();
  const int nhistobins = hSeed->GetNbinsX();
  float binmean[nfvbins][nhistobins+2];
  float binvar[nfvbins][nhistobins+2];
  
  // init to zero
  for (int fvbin=0; fvbin<nfvbins; fvbin++){
    for (int ibin=0; ibin<=nhistobins; ibin++){
      binmean[fvbin][ibin] = 0.;
      binvar[fvbin][ibin] = 0.;
    }
  }

  // add bin contents
  for (int fvbin=0; fvbin<nfvbins; fvbin++){
    // skip first throw (shold be default)
    for (int ithrow=1; ithrow<nPoints; ithrow++){
      for (int ibin=1; ibin<=nhistobins; ibin++){
        binmean[fvbin][ibin] += getHistogram(ithrow,fvbin)->GetBinContent(ibin);
      }
    }
  }

  // normalize
  float N=(float)(nPoints-1);
  for (int fvbin=0; fvbin<nfvbins; fvbin++){
    for (int ibin=0; ibin<=nhistobins; ibin++){
      binmean[fvbin][ibin] /= N;
    }
  }

  // fill histogram of systematic uncertainty in each bin of the
  // the 1D seed histograms
  for (int fvbin=0; fvbin<nfvbins; fvbin++){
    TString hname = "Bin_Uncertainty";
    hname.Append(Form("_FVBin%d",fvbin));
    binUnc[fvbin]=(TH1D*)hSeed->Clone(hname.Data());
    binUnc[fvbin]->Reset();
  }
  for (int fvbin=0; fvbin<nfvbins; fvbin++){
    for (int ibin=1; ibin<=nhistobins; ibin++){
      binUnc[fvbin]->SetBinContent(ibin,binmean[fvbin][ibin]);
    }
  }
 
  ////////////////////////////////////////////////////////
  // find variance in each bin
  for (int fvbin=0; fvbin<nfvbins; fvbin++){
    for (int ithrow=1; ithrow<nPoints; ithrow++){
      for (int ibin=1; ibin<=nhistobins; ibin++){
        binvar[fvbin][ibin] += ((getHistogram(ithrow,fvbin)->GetBinContent(ibin)-binmean[fvbin][ibin])*
                               (getHistogram(ithrow,fvbin)->GetBinContent(ibin)-binmean[fvbin][ibin]));
      }
    }
  }
  for (int fvbin=0; fvbin<nfvbins; fvbin++){
    for (int ibin=1; ibin<=nhistobins; ibin++){
      cout<<"var: "<<binvar[fvbin][ibin]<<endl;
      double error = TMath::Sqrt(binvar[fvbin][ibin]/N);
      binUnc[fvbin]->SetBinError(ibin,error);
    }
  }

  //////////////////////////////////////////////////////
  // total uncertainty maps in (towall,wall) space 
  FVUncMap = (TH2FV*)hFV[0]->Clone("FVUncMap");
  FVUncMapCCQE = (TH2FV*)hFV[0]->Clone("FVUncMapCCQE");
  FVUncMapCCnQE = (TH2FV*)hFV[0]->Clone("FVUncMapCCnQE");
  FVUncMapNC = (TH2FV*)hFV[0]->Clone("FVUncMapNC");
  FVUncMapCCWrong = (TH2FV*)hFV[0]->Clone("FVUncMapCCWrong");
  //
  FVShiftMap = (TH2FV*)hFV[0]->Clone("FVShiftMap");
  FVShiftMapCCQE = (TH2FV*)hFV[0]->Clone("FVShiftMapCCQE");
  FVShiftMapCCnQE = (TH2FV*)hFV[0]->Clone("FVShiftMapCCnQE");
  FVShiftMapCCWrong = (TH2FV*)hFV[0]->Clone("FVShiftMapCCWrong");
  FVShiftMapNC = (TH2FV*)hFV[0]->Clone("FVShiftMapNC");
  //
  FVFitMap = (TH2FV*)hFV[0]->Clone("FVFitMap");
  FVFitMapCCQE = (TH2FV*)hFV[0]->Clone("FVFitMapCCQE");
  FVFitMapCCnQE = (TH2FV*)hFV[0]->Clone("FVFitMapCCnQE");
  FVFitMapCCWrong = (TH2FV*)hFV[0]->Clone("FVFitMapCCWrong");
  FVFitMapNC = (TH2FV*)hFV[0]->Clone("FVFitMapNC");

  // calc RMS and mean for each event type in each FV bin
  float RMS[nfvbins][5]; 
  float NMean[nfvbins][5];
  for (int fvbin=0; fvbin<nfvbins; fvbin++){
    for (int ih=0; ih<5; ih++){
       RMS[fvbin][ih]=0.;
       NMean[fvbin][ih]=0.;
    }
  }
  // calculate mean
  float norm = (float)(nPoints-1);
  for (int fvbin=0; fvbin<nfvbins; fvbin++){
    for (int ithrow=1; ithrow<nPoints; ithrow++){
      //
      float nevts = hFV[ithrow]->GetBinContent(fvbin+1);
      NMean[fvbin][0] += (nevts)/norm;
      //
      nevts = hFVCCQE[ithrow]->GetBinContent(fvbin+1);
      NMean[fvbin][1] += (nevts)/norm;
      //
      nevts = hFVCCnQE[ithrow]->GetBinContent(fvbin+1);
      NMean[fvbin][2] += (nevts)/norm;
      //
      nevts = hFVCCWrong[ithrow]->GetBinContent(fvbin+1);
      NMean[fvbin][3] += (nevts)/norm;
      //
      nevts = hFVNC[ithrow]->GetBinContent(fvbin+1);
      NMean[fvbin][4] += (nevts)/norm;
      //
    }
  }
  // calculate variance
  for (int fvbin=0; fvbin<nfvbins; fvbin++){
    for (int ithrow=1; ithrow<nPoints; ithrow++){
      //
      float nevts = hFV[ithrow]->GetBinContent(fvbin+1);
      float n0 = NMean[fvbin][0];
      RMS[fvbin][0] += (nevts-n0)*(nevts-n0)/norm;
      //
      nevts = hFVCCQE[ithrow]->GetBinContent(fvbin+1);
      n0 = NMean[fvbin][1];
      RMS[fvbin][1] += (nevts-n0)*(nevts-n0)/norm;
      //
      nevts = hFVCCnQE[ithrow]->GetBinContent(fvbin+1);
      n0 = NMean[fvbin][2];
      RMS[fvbin][2] += (nevts-n0)*(nevts-n0)/norm;
      //
      nevts = hFVCCWrong[ithrow]->GetBinContent(fvbin+1);
      n0 = NMean[fvbin][3];
      RMS[fvbin][3] += (nevts-n0)*(nevts-n0)/norm;
      //
      nevts = hFVNC[ithrow]->GetBinContent(fvbin+1);
      n0 = NMean[fvbin][4];
      RMS[fvbin][4] += (nevts-n0)*(nevts-n0)/norm;
      //
    }
  }
  // make it RMS
  for (int fvbin=0; fvbin<nfvbins; fvbin++){
    for (int ih=0; ih<5; ih++){
       RMS[fvbin][ih]=TMath::Sqrt(RMS[fvbin][ih]);
    }
  }
 
  // make and fill histograms of the # of events in each FV bin
  // get max # of events from largest bin
  float Nmax = 0;;
  float width=0;
  float Nmin =0.;
  float factor = 5.4;
  int nbins = 80;
  // initial binning
  for (int fvbin=0; fvbin<nfvbins; fvbin++){
    //
    Nmax = NMean[fvbin][0]+(RMS[fvbin][0]*factor);
    Nmin = NMean[fvbin][0]-(RMS[fvbin][0]*factor);
    TString hname = "hNevents_";
    hname.Append(Form("%d",fvbin));
    hNevents[fvbin] = new TH1D(hname.Data(),hname.Data(),nbins,Nmin,Nmax);
    //
    Nmax = NMean[fvbin][1]+(RMS[fvbin][1]*factor);
    Nmin = NMean[fvbin][1]-(RMS[fvbin][1]*factor);
    hname = "hNeventsCCQE_";
    hname.Append(Form("%d",fvbin));
    hNeventsCCQE[fvbin] = new TH1D(hname.Data(),hname.Data(),nbins,Nmin,Nmax);
    //
    Nmax = NMean[fvbin][2]+(RMS[fvbin][2]*factor);
    Nmin = NMean[fvbin][2]-(RMS[fvbin][2]*factor);
    hname = "hNeventsCCnQE_";
    hname.Append(Form("%d",fvbin));
    hNeventsCCnQE[fvbin] = new TH1D(hname.Data(),hname.Data(),nbins,Nmin,Nmax);
    //
    Nmax = NMean[fvbin][3]+(RMS[fvbin][3]*factor);
    Nmin = NMean[fvbin][3]-(RMS[fvbin][3]*factor);
    hname = "hNeventsCCWrong_";
    hname.Append(Form("%d",fvbin));
    hNeventsCCWrong[fvbin] = new TH1D(hname.Data(),hname.Data(),nbins,Nmin,Nmax);
    //
    Nmax = NMean[fvbin][4]+(RMS[fvbin][4]*factor);
    Nmin = NMean[fvbin][4]-(RMS[fvbin][4]*factor);
    hname = "hNeventsNC";
    hname.Append(Form("%d",fvbin));
    hNeventsNC[fvbin] = new TH1D(hname.Data(),hname.Data(),nbins,Nmin,Nmax);
    //
  }
  // now fill values
  for (int fvbin=0; fvbin<nfvbins; fvbin++){
    for (int ithrow=1; ithrow<nPoints; ithrow++){
      //
      float nevts = hFV[ithrow]->GetBinContent(fvbin+1);
      hNevents[fvbin]->Fill(nevts);
      //
      nevts = hFVCCQE[ithrow]->GetBinContent(fvbin+1);
      hNeventsCCQE[fvbin]->Fill(nevts);
      //
      nevts = hFVCCnQE[ithrow]->GetBinContent(fvbin+1);
      hNeventsCCnQE[fvbin]->Fill(nevts);
      //
      nevts = hFVCCWrong[ithrow]->GetBinContent(fvbin+1);
      hNeventsCCWrong[fvbin]->Fill(nevts);
      //
      nevts = hFVNC[ithrow]->GetBinContent(fvbin+1);
      hNeventsNC[fvbin]->Fill(nevts);
      //
    }
  }

  //////////////////////////////////////////////////////////
  // fit to gaussians and calculate shifts
  fillFitAndShiftError();

  // fill lines at the values of the nominal contents
  nominalLine[0] = new TLine(hFV[0]->GetBinContent(1),0,hFV[0]->GetBinContent(1),10000);
  nominalLine[1] = new TLine(hFV[0]->GetBinContent(2),0,hFV[0]->GetBinContent(2),10000);
  nominalLine[2] = new TLine(hFV[0]->GetBinContent(3),0,hFV[0]->GetBinContent(3),1000);
  nominalLine[3] = new TLine(hFV[0]->GetBinContent(4),0,hFV[0]->GetBinContent(4),10000);
  nominalLine[4] = new TLine(hFV[0]->GetBinContent(5),0,hFV[0]->GetBinContent(5),10000);
  nominalLine[5] = new TLine(hFV[0]->GetBinContent(6),0,hFV[0]->GetBinContent(6),10000);

  // find the total uncertainty
  for (int fvbin=0; fvbin<hFV[0]->GetNumberOfBins(); fvbin++){
    //
    float binc1 = FVFitMap->GetBinContent(fvbin+1); //< this is "fit" error
    float binc2 = FVShiftMap->GetBinContent(fvbin+1); //< this is "shift" error
    FVUncMap->SetBinContent(fvbin+1,TMath::Sqrt(binc1*binc1 + binc2*binc2));
    //
    binc1 = FVFitMapCCQE->GetBinContent(fvbin+1);
    binc2 = FVShiftMapCCQE->GetBinContent(fvbin+1);
    FVUncMapCCQE->SetBinContent(fvbin+1,TMath::Sqrt(binc1*binc1 + binc2*binc2));
    //
    binc1 = FVFitMapCCnQE->GetBinContent(fvbin+1);
    binc2 = FVShiftMapCCnQE->GetBinContent(fvbin+1);
    FVUncMapCCnQE->SetBinContent(fvbin+1,TMath::Sqrt(binc1*binc1 + binc2*binc2));
    //
    binc1 = FVFitMapCCWrong->GetBinContent(fvbin+1);
    binc2 = FVShiftMapCCWrong->GetBinContent(fvbin+1);
    FVUncMapCCWrong->SetBinContent(fvbin+1,TMath::Sqrt(binc1*binc1 + binc2*binc2));
    //
    binc1 = FVFitMapNC->GetBinContent(fvbin+1);
    binc2 = FVShiftMapNC->GetBinContent(fvbin+1);
    FVUncMapNC->SetBinContent(fvbin+1,TMath::Sqrt(binc1*binc1 + binc2*binc2));
  }

  //
  return;
 
}

////////////////////////////////////////////////////////////
// fill Fit and Shift Erros
void modHistoArrayFV::fillFitAndShiftError(){

  // make array of gaussians
  const int nfvbins = hFV[0]->GetNumberOfBins();
  // for total # of events
  for (int fvbin=1; fvbin<=nfvbins; fvbin++){

    // Total
    float Norm =  hFV[0]->GetBinContent(fvbin);
    float fractional_uncertainty = hNevents[fvbin-1]->GetRMS()/Norm;
    float fractional_shift = ( hNevents[fvbin-1]->GetMean()- hFV[0]->GetBinContent(fvbin));
    fractional_shift /= Norm;
    FVFitMap->SetBinContent(fvbin,fractional_uncertainty);
    FVShiftMap->SetBinContent(fvbin,fractional_shift);
    // CCQE
    Norm =  hFVCCQE[0]->GetBinContent(fvbin);
    fractional_uncertainty = hNeventsCCQE[fvbin-1]->GetRMS()/Norm;
    fractional_shift = (hNeventsCCQE[fvbin-1]->GetMean() - hFVCCQE[0]->GetBinContent(fvbin));
    fractional_shift /= Norm;
    FVFitMapCCQE->SetBinContent(fvbin,fractional_uncertainty);
    FVShiftMapCCQE->SetBinContent(fvbin,fractional_shift);    
    // CCnQE
    Norm =  hFVCCnQE[0]->GetBinContent(fvbin);
    fractional_uncertainty = hNeventsCCnQE[fvbin-1]->GetRMS()/Norm;
    fractional_shift = (hNeventsCCnQE[fvbin-1]->GetMean() - hFVCCnQE[0]->GetBinContent(fvbin));
    fractional_shift /= Norm;
    FVFitMapCCnQE->SetBinContent(fvbin,fractional_uncertainty);
    FVShiftMapCCnQE->SetBinContent(fvbin,fractional_shift);
    // CCWrong
    Norm = hFVCCWrong[0]->GetBinContent(fvbin);
    fractional_uncertainty = hNeventsCCWrong[fvbin-1]->GetRMS()/Norm;
    fractional_shift = (hNeventsCCWrong[fvbin-1]->GetMean() - hFVCCWrong[0]->GetBinContent(fvbin));
    fractional_shift /= Norm;
    FVFitMapCCWrong->SetBinContent(fvbin,fractional_uncertainty);
    FVShiftMapCCWrong->SetBinContent(fvbin,fractional_shift);
    // NC
    Norm =  hFVNC[0]->GetBinContent(fvbin);
    fractional_uncertainty = hNeventsNC[fvbin-1]->GetRMS()/Norm;
    fractional_shift = (hNeventsNC[fvbin-1]->GetMean() - hFVNC[0]->GetBinContent(fvbin));
    fractional_shift /= Norm;
    FVFitMapNC->SetBinContent(fvbin,fractional_uncertainty);
    FVShiftMapNC->SetBinContent(fvbin,fractional_shift);

  }


  
  return; 
}



////////////////////////////////////////////////////////////
// print to directory
void modHistoArrayFV::printUncMap(const char* plotdir){

  TCanvas* cc = new TCanvas("cc","cc",700,800);
  cc->Divide(2,3);
 
  // number of events
  for (int fvbin=0; fvbin<hFV[0]->GetNumberOfBins(); fvbin++){
    cc->cd(fvbin+1);
    double mean = hNevents[fvbin]->GetMean();
    double std_dev = hNevents[fvbin]->GetRMS();
    double xmin = mean-(6*std_dev);
    double xmax = mean+(6*std_dev);
    double nominal = hFV[0]->GetBinContent(fvbin+1);
    nominalLine[fvbin] = new TLine(nominal,0,nominal,1000);
    hNevents[fvbin]->GetXaxis()->SetRangeUser(xmin,xmax);
    hNevents[fvbin]->SetBit(TH1::kNoTitle);
    hNevents[fvbin]->GetXaxis()->SetTitle("# of events");
    hNevents[fvbin]->GetXaxis()->SetTitleSize(0.05);
    hNevents[fvbin]->GetYaxis()->SetTitle("# of throws");
    hNevents[fvbin]->GetYaxis()->SetTitleSize(0.05);
    hNevents[fvbin]->Draw();
  }
  TString plotname = plotdir;
  plotname.Append("NumOfEvents_Unc.png");
  cc->Print(plotname.Data());

  // number of events ccqe
  for (int fvbin=0; fvbin<hFV[0]->GetNumberOfBins(); fvbin++){
    cc->cd(fvbin+1);
    double mean = hNeventsCCQE[fvbin]->GetMean();
    double std_dev = hNeventsCCQE[fvbin]->GetRMS();
    double xmin = mean-(7*std_dev);
    double xmax = mean+(7*std_dev);
    double nominal = hFVCCQE[0]->GetBinContent(fvbin+1);
    nominalLine[fvbin] = new TLine(nominal,0,nominal,1000);
    hNeventsCCQE[fvbin]->GetXaxis()->SetRangeUser(xmin,xmax);
    hNeventsCCQE[fvbin]->SetBit(TH1::kNoTitle);
    hNeventsCCQE[fvbin]->GetXaxis()->SetTitle("# of events");
    hNeventsCCQE[fvbin]->GetXaxis()->SetTitleSize(0.05);
    hNeventsCCQE[fvbin]->GetYaxis()->SetTitle("# of throws");
    hNeventsCCQE[fvbin]->GetYaxis()->SetTitleSize(0.05);
    hNeventsCCQE[fvbin]->Draw();
    nominalLine[fvbin]->SetLineColor(kBlue);
    nominalLine[fvbin]->Draw("same");
  }
  plotname = plotdir;
  plotname.Append("NumOfEventsCCQE_Unc.png");
  cc->Print(plotname.Data());

  // number of events ccnqe
  for (int fvbin=0; fvbin<hFV[0]->GetNumberOfBins(); fvbin++){
    cc->cd(fvbin+1);
    double mean = hNeventsCCnQE[fvbin]->GetMean();
    double std_dev = hNeventsCCnQE[fvbin]->GetRMS();
    double xmin = mean-(7*std_dev);
    double xmax = mean+(7*std_dev);
    double nominal = hFVCCnQE[0]->GetBinContent(fvbin+1);
    nominalLine[fvbin] = new TLine(nominal,0,nominal,1000);
    hNeventsCCnQE[fvbin]->GetXaxis()->SetRangeUser(xmin,xmax);
    hNeventsCCnQE[fvbin]->SetBit(TH1::kNoTitle);
    hNeventsCCnQE[fvbin]->GetXaxis()->SetTitle("# of events");
    hNeventsCCnQE[fvbin]->GetXaxis()->SetTitleSize(0.05);
    hNeventsCCnQE[fvbin]->GetYaxis()->SetTitle("# of throws");
    hNeventsCCnQE[fvbin]->GetYaxis()->SetTitleSize(0.05);
    hNeventsCCnQE[fvbin]->Draw();
    nominalLine[fvbin]->SetLineColor(kBlue);
    nominalLine[fvbin]->Draw("same");
  }
  plotname = plotdir;
  plotname.Append("NumOfEventsCCnQE_Unc.png");
  cc->Print(plotname.Data());


  // number of events ccwrong 
  for (int fvbin=0; fvbin<hFV[0]->GetNumberOfBins(); fvbin++){
    cc->cd(fvbin+1);
    double mean = hNeventsCCWrong[fvbin]->GetMean();
    double std_dev = hNeventsCCWrong[fvbin]->GetRMS();
    double xmin = mean-(7*std_dev);
    double xmax = mean+(7*std_dev);
    double nominal = hFVCCWrong[0]->GetBinContent(fvbin+1);
    nominalLine[fvbin] = new TLine(nominal,0,nominal,1000);
    hNeventsCCWrong[fvbin]->GetXaxis()->SetRangeUser(xmin,xmax);
    hNeventsCCWrong[fvbin]->SetBit(TH1::kNoTitle);
    hNeventsCCWrong[fvbin]->GetXaxis()->SetTitle("# of events");
    hNeventsCCWrong[fvbin]->GetXaxis()->SetTitleSize(0.05);
    hNeventsCCWrong[fvbin]->GetYaxis()->SetTitle("# of throws");
    hNeventsCCWrong[fvbin]->GetYaxis()->SetTitleSize(0.05);
    hNeventsCCWrong[fvbin]->Draw();
    nominalLine[fvbin]->SetLineColor(kBlue);
    nominalLine[fvbin]->Draw("same");
  }
  plotname = plotdir;
  plotname.Append("NumOfEventsCCWrong_Unc.png");
  cc->Print(plotname.Data());



  // number of events nc 
  for (int fvbin=0; fvbin<hFV[0]->GetNumberOfBins(); fvbin++){
    cc->cd(fvbin+1);
    double mean = hNeventsNC[fvbin]->GetMean();
    double std_dev = hNeventsNC[fvbin]->GetRMS();
    double xmin = mean-(7*std_dev);
    double xmax = mean+(7*std_dev);
    double nominal = hFVNC[0]->GetBinContent(fvbin+1);
    nominalLine[fvbin] = new TLine(nominal,0,nominal,1000);
    hNeventsNC[fvbin]->GetXaxis()->SetRangeUser(xmin,xmax);
    hNeventsNC[fvbin]->SetBit(TH1::kNoTitle);
    hNeventsNC[fvbin]->GetXaxis()->SetTitle("# of events");
    hNeventsNC[fvbin]->GetXaxis()->SetTitleSize(0.05);
    hNeventsNC[fvbin]->GetYaxis()->SetTitle("# of throws");
    hNeventsNC[fvbin]->GetYaxis()->SetTitleSize(0.05);
    hNeventsNC[fvbin]->Draw();
    nominalLine[fvbin]->SetLineColor(kBlue);
    nominalLine[fvbin]->Draw("same");
  }
  plotname = plotdir;
  plotname.Append("NumOfEventsNC_Unc.png");
  cc->Print(plotname.Data());

  // distribution uncertainty
  for (int fvbin=0; fvbin<hFV[0]->GetNumberOfBins(); fvbin++){
    cc->cd(fvbin+1);
    binUnc[fvbin]->SetBit(TH1::kNoTitle);    
    hNevents[fvbin]->SetStats(0);
    binUnc[fvbin]->SetFillColor(kOrange);
    binUnc[fvbin]->GetXaxis()->SetTitle("Erec [MeV]");
    binUnc[fvbin]->GetYaxis()->SetRangeUser(0.,binUnc[fvbin]->GetMaximum()*2);
    binUnc[fvbin]->Draw("e2");
    histos[0][fvbin]->SetLineColor(kBlue);
    histos[0][fvbin]->SetLineWidth(3);
    histos[0][fvbin]->Draw("same");
    nominalLine[fvbin]->SetLineColor(kBlue);
    nominalLine[fvbin]->Draw("same");
  }
  plotname = plotdir;
  plotname.Append("Enu_Distributions.png");
  cc->Print(plotname.Data());


}


///////////////////////////////////////////////////////////
// constructor from seed histogram, 
modHistoArrayFV::modHistoArrayFV(TH1D* hseed, TH2FV* hfv, int ninit){

  // get base name tag
  cout<<"modHistoArrayFV: Called constructor"<<nameTag.Data()<<endl;
  nameTag = hseed->GetName();
  cout<<"modHistoArrayFV: Constructing with name"<<nameTag.Data()<<endl;
 
  // output file
  TString fname = nameTag.Data();
  fname.Append("_histogram_array.root");
  cout<<"modHistoArrayFV: Opening file: "<<nameTag.Data()<<endl;
  fout = new TFile(fname.Data(),"RECREATE");

  // setup histogram seeds
  TString hname = nameTag.Data();
  hname.Append("_seed");
  hSeed = (TH1D*)hseed->Clone(hname.Data());
  TString hfvname = nameTag.Data();
  hfvname.Append("_fvseed");
  hFV[0] = (TH2FV*)hfv->Clone(hfvname.Data());
  
  // initialize array
  nHistos = 0;
  nPoints = 0;
  currentIndex = 0;
  if (ninit>=0) init(ninit);

}

///////////////////////////////////////////////////////////////////////////////////
// read from file
void modHistoArrayFV::readFromFile(int nfvbins, int nthrows){
 
  TList* keylist = fout->GetListOfKeys();

  // for histos[] array
  for (int fvbin=0; fvbin<nfvbins; fvbin++){
    for (int ithrow=0; ithrow<nthrows; ithrow++){
      TString wantkey = Form("throw%d_fvbin%d",ithrow,fvbin);
      cout<<"Looking for key with: "<<wantkey.Data()<<endl;
      for (int ikey=0; ikey<keylist->GetSize(); ikey++){
        TString keycheck = keylist->At(ikey)->GetName();
        if (keycheck.Contains(wantkey.Data())){
          cout<<"Found "<<keycheck.Data()<<"!"<<endl;
          histos[ithrow][fvbin] = (TH1D*)fout->Get(keycheck.Data());
          break;
        }
      }
    }
  }
 

  // get seed array
  TString wantkey = "_seed";  
  for (int ikey=0; ikey<keylist->GetSize(); ikey++){
    cout<<"Looking for key with: "<<wantkey.Data()<<endl;
    TString keycheck = keylist->At(ikey)->GetName();
    if (keycheck.Contains(wantkey.Data())){
      cout<<"Found "<<keycheck.Data()<<"!"<<endl;
      hSeed = (TH1D*)fout->Get(keycheck.Data());
      break;
    }
  }  

  // for hFV[] array
  wantkey = "fvseed";  
  for (int ikey=0; ikey<keylist->GetSize(); ikey++){
    cout<<"Looking for key with: "<<wantkey.Data()<<endl;
    TString keycheck = keylist->At(ikey)->GetName();
    if (keycheck.Contains(wantkey.Data())){
      cout<<"Found "<<keycheck.Data()<<"!"<<endl;
      hFV[0] = (TH2FV*)fout->Get(keycheck.Data());
      break;
    }
  }  

  for (int ithrow=0; ithrow<nthrows; ithrow++){
    wantkey = Form("FV_throw%d",ithrow);
    for (int ikey=0; ikey<keylist->GetSize(); ikey++){
      TString keycheck = keylist->At(ikey)->GetName();
      if (keycheck.Contains(wantkey.Data())){
        if (keycheck.Contains("CCQE")){
          cout<<"Found "<<keycheck.Data()<<"!"<<endl;
          hFVCCQE[ithrow] = (TH2FV*)fout->Get(keycheck.Data());
//          break;
        }
        else if (keycheck.Contains("CCnQE")){
          cout<<"Found "<<keycheck.Data()<<"!"<<endl;
          hFVCCnQE[ithrow] = (TH2FV*)fout->Get(keycheck.Data());
//          break;
        }
        else if (keycheck.Contains("CCWrong")){
          cout<<"Found "<<keycheck.Data()<<"!"<<endl;
          hFVCCWrong[ithrow] = (TH2FV*)fout->Get(keycheck.Data());
//          break;
        }
        else if (keycheck.Contains("NC")){
          cout<<"Found "<<keycheck.Data()<<"!"<<endl;
          hFVNC[ithrow] = (TH2FV*)fout->Get(keycheck.Data());
//          break;
        }
        else{
          cout<<"Found "<<keycheck.Data()<<"!"<<endl;
          hFV[ithrow] = (TH2FV*)fout->Get(keycheck.Data());
        }
//        break;
      }
    }
  }
  

  return;

}

///////////////////////////////////////////////////////////////////////////////////
// Construct from file
modHistoArrayFV::modHistoArrayFV(const char* filename, int nfvbins, int nthrows){
  
  fout = new TFile(filename);
  nPoints = nthrows;
  nHistos = 0;

  // read in histograms frm file
  readFromFile(nfvbins,nthrows);

}


///////////////////////////////////
// get a specific histogram
TH1D* modHistoArrayFV::getHistogram(int ithrow, int fvbin){
  return histos[ithrow][fvbin];
}

//////////////////////////////////////
// draw the distribution of Nev for a given bin
void modHistoArrayFV::drawNev(int fvbin){
  hNevents[fvbin]->Draw();
  nominalLine[fvbin]->SetLineColor(kBlue);
  nominalLine[fvbin]->Draw("same");

  
}

///////////////////////////////////
//draw all histos
void modHistoArrayFV::drawArray(int fvbin){
  
  histos[0][fvbin]->SetLineColor(kBlue);
  histos[0][fvbin]->Draw();
  for (int i=0; i<nPoints; i++){
    histos[i][fvbin]->Draw("same");
  }

  return;
}

///////////////////////////////////
//setup array of histos
void modHistoArrayFV::init(int ninit){

  nHistos = 0;
  nPoints = ninit;

  // initialize all histos
  // (clone them from hFV[0])
  for (int i=0; i<ninit; i++){

       // all histos
       TString hnamebase = nameTag.Data();
       hnamebase.Append(Form("_FV_throw%d",i)); 
       if (i>0) hFV[i] = (TH2FV*)hFV[0]->Clone(hnamebase.Data());

       // ccqe
       hnamebase = nameTag.Data();
       hnamebase.Append("_CCQE");
       hnamebase.Append(Form("_FV_throw%d",i)); 
       hFVCCQE[i] = (TH2FV*)hFV[0]->Clone(hnamebase.Data());

       // ccnqe
       hnamebase = nameTag.Data();
       hnamebase.Append("_CCnQE");
       hnamebase.Append(Form("_FV_throw%d",i)); 
       hFVCCnQE[i] = (TH2FV*)hFV[0]->Clone(hnamebase.Data());

       // ccwrong flavor
       hnamebase = nameTag.Data();
       hnamebase.Append("_CCWrong");
       hnamebase.Append(Form("_FV_throw%d",i)); 
       hFVCCWrong[i] = (TH2FV*)hFV[0]->Clone(hnamebase.Data());

       // nc
       hnamebase = nameTag.Data();
       hnamebase.Append("_NC");
       hnamebase.Append(Form("_FV_throw%d",i)); 
       hFVNC[i] = (TH2FV*)hFV[0]->Clone(hnamebase.Data());

    // initialize all 1D histos
    for (int fvbin=0; fvbin<hFV[0]->GetNumberOfBins(); fvbin++){
      TString hname = nameTag.Data();
      hname.Append(Form("_throw%d",i));
      hname.Append(Form("_fvbin%d",fvbin));
      histos[i][fvbin] = (TH1D*)hSeed->Clone(hname.Data());
      histos[i][fvbin]->Reset();
      nHistos++;
    }
  } 

  currentIndex = 0;

  return;
}

///////////////////////////////////////////////////
//copy hisogram contents to a point in the array
void modHistoArrayFV::setHistoContents(TH1D* hadd, int index, int fvbin){

//  if (index<0){
//    index = currentIndex;
//    nHistos++;
//  }
 
  int nbins = histos[index][fvbin]->GetNbinsX();
  for (int ibin=1; ibin<=nbins; ibin++){
    double binc = hadd->GetBinContent(ibin);
    double binerr = hadd->GetBinError(ibin);
    histos[index][fvbin]->SetBinContent(ibin,binc);
    histos[index][fvbin]->SetBinError(ibin,binerr);
  }
 
  currentIndex++;
  
  return;
}



#endif





