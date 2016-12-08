#ifndef HISOTARRAYFV_H
#define HISOTARRAYFV_H

#include "modHistoArrayFV.h"



using namespace std;

////////////////////////////////////////////////////////////////////////
// Save the overall uncertainties to be used in moreUncertainties.cxx
void modHistoArrayFV::saveSummary(const char* dir){

  // set up file
  TString outfilename = dir;
  outfilename.Append("FVUncMap.root");
  TFile* outfile = new TFile(outfilename.Data(),"RECREATE");
  outfile->cd();

  // write histograms
  FVUncMap->Write();
  FVShiftMap->Write();

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

  outfile->Close();
  
  //
  return;
}



///////////////////////////////////////
// save histograms and quit
void modHistoArrayFV::saveClose(){
  fout->Write();
  fout->Close();
  return;
}



///////////////////////////////////////////////////////////
// Calculate summary statistics
void modHistoArrayFV::calcSummary(){

  /////////////////////////////////////////////////////
  // find mean for each bin
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
  // total uncertainty map 
  FVUncMap = (TH2FV*)hFV[0]->Clone("FVUncMap");
  FVShiftMap = (TH2FV*)hFV[0]->Clone("FVShiftMap");
  float Nmax = hFV[0]->GetBinContent(nfvbins)*1.3;
//  float Nmin = hFV[0]->GetBinContent(nfvbins)*1.3;
  int nbins = 1900;
  for (int fvbin=0; fvbin<nfvbins; fvbin++){
    TString hname = "hNevents_";
    hname.Append(Form("%d",fvbin));
//    float Nmax = hFV[0]->GetBinContent(fvbin+1)*(1.5);
//    float Nmin = hFV[0]->GetBinContent(fvbin+1)*(0.5);
    hNevents[fvbin] = new TH1D(hname.Data(),hname.Data(),nbins,0,Nmax);
  }
  for (int fvbin=0; fvbin<nfvbins; fvbin++){
    for (int ithrow=1; ithrow<nPoints; ithrow++){
      float nevts = hFV[ithrow]->GetBinContent(fvbin+1);
      hNevents[fvbin]->Fill(nevts);
    }
  }

  //////////////////////////////////////////////////////////
  // fit to gaussians and calculate shifts
  fitGaussians();

  // fill lines at the values of the nominal contents
  nominalLine[0] = new TLine(hFV[0]->GetBinContent(1),0,hFV[0]->GetBinContent(1),10000);
  nominalLine[1] = new TLine(hFV[0]->GetBinContent(2),0,hFV[0]->GetBinContent(2),10000);
  nominalLine[2] = new TLine(hFV[0]->GetBinContent(3),0,hFV[0]->GetBinContent(3),1000);
  nominalLine[3] = new TLine(hFV[0]->GetBinContent(4),0,hFV[0]->GetBinContent(4),10000);
  nominalLine[4] = new TLine(hFV[0]->GetBinContent(5),0,hFV[0]->GetBinContent(5),10000);
  nominalLine[5] = new TLine(hFV[0]->GetBinContent(6),0,hFV[0]->GetBinContent(6),10000);

  // find the total uncertainty
  for (int fvbin=0; fvbin<hFV[0]->GetNumberOfBins(); fvbin++){
    float binc1 = FVUncMap->GetBinContent(fvbin+1);
    float binc2 = FVShiftMap->GetBinContent(fvbin+1);
    FVUncMap->SetBinContent(fvbin+1,binc1+binc2);
  }


  return;
 
}

////////////////////////////////////////////////////////////
// Fit Gaussians to Nevent distributions
void modHistoArrayFV::fitGaussians(){

  // make array of gaussians
  const int nfvbins = hFV[0]->GetNumberOfBins();
  TF1* gaussians[nfvbins];
  for (int fvbin=0; fvbin<nfvbins; fvbin++){
    TString fname = Form("fgauss_%d",fvbin);
    gaussians[fvbin] = new TF1(fname.Data(),"gaus",0,
                               hNevents[0]->GetBinLowEdge(hNevents[0]->GetNbinsX())*1.5);
    gaussians[fvbin]->SetParameter(0,hNevents[fvbin]->GetMaximum());
    gaussians[fvbin]->SetParameter(1,hNevents[fvbin]->GetMean());
    gaussians[fvbin]->SetParameter(2,hNevents[fvbin]->GetRMS());
    gaussians[fvbin]->SetLineColor(kRed);
    hNevents[fvbin]->Fit(fname.Data());
    float fractional_uncertainty = (gaussians[fvbin]->GetParameter(2)/gaussians[fvbin]->GetParameter(1));
    float fractional_shift = TMath::Abs(gaussians[fvbin]->GetParameter(1)- hFV[0]->GetBinContent(fvbin+1));
    fractional_shift /= gaussians[fvbin]->GetParameter(1);
    FVUncMap->SetBinContent(fvbin+1,fractional_uncertainty);
    FVShiftMap->SetBinContent(fvbin+1,fractional_shift);
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
    hNevents[fvbin]->GetXaxis()->SetRangeUser(xmin,xmax);
    hNevents[fvbin]->SetBit(TH1::kNoTitle);
    hNevents[fvbin]->GetXaxis()->SetTitle("# of events");
    hNevents[fvbin]->GetXaxis()->SetTitleSize(0.05);
    hNevents[fvbin]->GetYaxis()->SetTitle("# of throws");
    hNevents[fvbin]->GetYaxis()->SetTitleSize(0.05);
    hNevents[fvbin]->Draw();
    TString plotname = plotdir;
    plotname.Append("NumOfEvents_Unc.png");
    cc->Print(plotname.Data());
  }
  TString plotname = plotdir;
  plotname.Append("NumOfEvents_Unc.png");
  cc->Print(plotname.Data());

  // distribution uncertainty
  for (int fvbin=0; fvbin<hFV[0]->GetNumberOfBins(); fvbin++){
    cc->cd(fvbin+1);
    binUnc[fvbin]->SetBit(TH1::kNoTitle);    
    hNevents[fvbin]->SetStats(0);
    binUnc[fvbin]->SetFillColor(kOrange);
    binUnc[fvbin]->GetXaxis()->SetTitle("Erec [MeV]");
    binUnc[fvbin]->Draw("e2");
    histos[0][fvbin]->SetLineColor(kBlue);
    histos[0][fvbin]->Draw("same");
  }
  plotname = plotdir;
  plotname.Append("Enu_Distributions.png");
  cc->Print(plotname.Data());


}


///////////////////////////////////////////////////////////
// constructor from seed histogram, 
modHistoArrayFV::modHistoArrayFV(TH1D* hseed, TH2FV* hfv, int ninit){

  // get base name tag
  nameTag = hseed->GetName();
 
  // output file
  TString fname = nameTag.Data();
  fname.Append("_histogram_array.root");
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

  for (int ithrow=1; ithrow<nthrows; ithrow++){
    wantkey = Form("FV_throw%d",ithrow);
    for (int ikey=0; ikey<keylist->GetSize(); ikey++){
      cout<<"Looking for key with: "<<wantkey.Data()<<endl;
      TString keycheck = keylist->At(ikey)->GetName();
      if (keycheck.Contains(wantkey.Data())){
        cout<<"Found "<<keycheck.Data()<<"!"<<endl;
        hFV[ithrow] = (TH2FV*)fout->Get(keycheck.Data());
        break;
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
  for (int i=0; i<ninit; i++){
    if (i>0){
       TString hname = nameTag.Data();
       hname.Append(Form("_FV_throw%d",i)); 
       hFV[i] = (TH2FV*)hFV[0]->Clone(hname.Data());
    }
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





