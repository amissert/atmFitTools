#ifndef HISOTARRAYFV_H
#define HISOTARRAYFV_H

#include "modHistoArrayFV.h"



using namespace std;


void modHistoArrayFV::save(){
  fout->Write();
  return;
}




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
    for (int ithrow=0; ithrow<nPoints; ithrow++){
      for (int ibin=1; ibin<=nhistobins; ibin++){
        binmean[fvbin][ibin] += getHistogram(ithrow,fvbin)->GetBinContent(ibin);
      }
    }
  }

  // normalize
  float N=(float)nPoints;
  for (int fvbin=0; fvbin<nfvbins; fvbin++){
    for (int ibin=0; ibin<=nhistobins; ibin++){
      binmean[fvbin][ibin] /= N;
    }
  }
  // fill summary histogram
  for (int fvbin=0; fvbin<nfvbins; fvbin++){
    TString hname = "bin_uncertainty";
    hname.Append(Form("_fvbin%d",fvbin));
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
    for (int ithrow=0; ithrow<nPoints; ithrow++){
      for (int ibin=1; ibin<=nhistobins; ibin++){
        binvar[fvbin][ibin] += ((getHistogram(ithrow,fvbin)->GetBinContent(ibin)-binmean[fvbin][ibin])*
                               (getHistogram(ithrow,fvbin)->GetBinContent(ibin)-binmean[fvbin][ibin]));
      }
    }
  }
  for (int fvbin=0; fvbin<nfvbins; fvbin++){
    for (int ibin=1; ibin<=nhistobins; ibin++){
      cout<<"var: "<<binvar[fvbin][ibin]<<endl;
      binUnc[fvbin]->SetBinError(ibin,TMath::Sqrt(binvar[fvbin][ibin]/N));
    }
  }

  //////////////////////////////////////////////////////
  // total uncertainty map 
  FVUncMap = (TH2FV*)hFV[0]->Clone("FVUncMap");
  float Nmax = hFV[0]->GetBinContent(nfvbins)*1.3;
  for (int fvbin=0; fvbin<nfvbins; fvbin++){
    TString hname = "hNevents_";
    hname.Append(Form("%d",fvbin));
    hNevents[fvbin] = new TH1D(hname.Data(),hname.Data(),800,0,Nmax);
  }
  for (int fvbin=0; fvbin<nfvbins; fvbin++){
    for (int ithrow=0; ithrow<nPoints; ithrow++){
      float nevts = hFV[ithrow]->GetBinContent(fvbin+1);
      hNevents[fvbin]->Fill(nevts);
    }
  }

  //////////////////////////////////////////////////////////
  // fit to gaussians
  fitGaussians();

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
    float fractional_uncertainty = 100.*(gaussians[fvbin]->GetParameter(2)/gaussians[fvbin]->GetParameter(1));
    FVUncMap->SetBinContent(fvbin+1,fractional_uncertainty);
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
  }
  plotname = plotdir;
  plotname.Append("Enu_Distributions.png");
  cc->Print(plotname.Data());


}


///////////////////////////////////////////////////////////
// constructor
modHistoArrayFV::modHistoArrayFV(TH1D* hseed, TH2FV* hfv, int ninit){

  // get base name tag
  nameTag = hseed->GetName();
 
  // output file
  TString fname = nameTag.Data();
  fname.Append("_histos.root");
  fout = new TFile(fname.Data(),"RECREATE");

  // setup histogram seed
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

///////////////////////////////////
//draw all histos
void modHistoArrayFV::drawArray(int fvbin){
  
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





