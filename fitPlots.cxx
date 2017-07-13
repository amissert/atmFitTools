#ifndef FITPLOTS_C
#define FITPLOTS_C

#include "fitPlots.h"
#include "hBound.cxx"



/////////////////////////////////////////////////////////////
void fitPlots::printFitSummaryFV(const char* dir, const char* name){

  // setup names
  TString nametag = name;
  TString outdir = dir;
  TString plotname = "";
  
  // total number of samples and bins
  int natt  = 4;
  int nsamp = 3;

  // make and print plot for each sample and bin
  for (int iatt=0; iatt<natt; iatt++){
    for (int isamp=0; isamp<nsamp; isamp++){
     
     // plot name
     plotname = outdir.Data();
     plotname.Append(nametag.Data());
     plotname.Append(Form("_samp%d_attribute%d.png",isamp,iatt));

     // draw
     drawFitSummaryFV(isamp, iatt);

     // print
     cc->Print(plotname.Data());
     cc->Delete();
    }
  }

  //
  return;
}


/////////////////////////////////////////////////////////////
void fitPlots::printFitSummary(const char* dir, const char* name){

  // setup names
  TString nametag = name;
  TString outdir = dir;
  TString plotname = "";
  
  // total number of samples and bins
  int nbins = 6;
  int nsamp = 3;

  // make and print plot for each sample and bin
  for (int ibin=0; ibin<nbins; ibin++){
    for (int isamp=0; isamp<nsamp; isamp++){
     
     // plot name
     plotname = outdir.Data();
     plotname.Append(nametag.Data());
     plotname.Append(Form("_samp%d_bin%d.png",isamp,ibin));

     // draw
     drawFitSummary(isamp, ibin);

     // print
     cc->Print(plotname.Data());
    
    }
  }

  //
  return;
}



/////////////////////////////////////////////////////////////
fitPlots::fitPlots(histoCompare* hc, TTree* partree){
   
   nAtt = 4;
   nPoints = 10;
   hCompare = hc;
   parTree = partree;
   mcmcPars = new mcmcReader(parTree);
   
   
}



/////////////////////////////////////////////////////////////
TString fitPlots::getAxisTitle(int iatt){

  // setup titles
  const int ntitles = 4;
  if (iatt > ntitles){
    cout<<"fitPlots::getAxisTitle: no title found!"<<endl;
    TString title = "none";
    return title;
  }
  TString xtitle[] = {"e/#mu PID",
                      "e/#pi^{0} PID",
                      "#mu/#pi PID",
                      "RC Parameter"};

  TString title = xtitle[iatt];
  return title;

}



/////////////////////////////////////////////////////////////
void fitPlots::drawFitSummaryFV(int isamp , int iatt){


  // setup titles
  TString xtitle[] = {"e/#mu PID",
                      "e/#pi^{0} PID",
                      "#mu/#pi PID",
                      "RC Parameter"};
  TString ytitle = "# of events";
  double titlesize = 0.05;
  double yoffset = 1.3;
  int axisdiv = 10;
  double binbuffer = 3.;

  cc = new TCanvas("cc","cc",1200,800);
  cc->Divide(3,2);

  // fill arrays
  const int nfvbins = 6;
  for (int ibin=0; ibin<nfvbins; ibin++){

    // run toy experiments
    fillArrays(ibin, isamp);

    // overall title
    TString title = Form("Detector Region %d",ibin);

    // go to proper pad
    int ipad = ibin+1;
    cc->cd(ipad);
    cc->GetPad(ipad)->SetLeftMargin(0.15);
    cc->GetPad(ipad)->SetRightMargin(0.01);
    cc->GetPad(ipad)->SetTopMargin(0.10);
    cc->GetPad(ipad)->SetBottomMargin(0.10);

    // setup mcmc summary
    hAtt[iatt]->calcSummary();
    hAtt[iatt]->hSummary->SetMarkerStyle(4);
    hAtt[iatt]->hSummary->GetXaxis()->SetTitle(getAxisTitle(iatt).Data());
    hAtt[iatt]->hSummary->GetYaxis()->SetTitle(ytitle.Data());
    hAtt[iatt]->hSummary->GetXaxis()->SetTitleSize(titlesize);
    hAtt[iatt]->hSummary->GetXaxis()->SetNdivisions(axisdiv);
    hAtt[iatt]->hSummary->GetYaxis()->SetNdivisions(axisdiv);
    hAtt[iatt]->hSummary->SetTitle(title.Data());

    // setup MC original
    TH1D* hmc = hCompare->hManager->getSumHistogram(isamp,ibin,iatt,1);
    hmc->SetLineColor(kRed);

    // setup data
    TH1D* hd = hCompare->hManager->hData[isamp][ibin][iatt];
    hd->SetMarkerStyle(8);
    hd->GetXaxis()->SetTitle(getAxisTitle(iatt).Data());
    hd->GetYaxis()->SetTitle(ytitle.Data());
    hd->GetYaxis()->SetTitleOffset(yoffset);
    hd->GetYaxis()->SetTitleSize(titlesize);
    hd->GetXaxis()->SetTitleSize(titlesize);
    hd->GetXaxis()->SetNdivisions(axisdiv);
    hd->GetYaxis()->SetNdivisions(axisdiv);
    hd->SetTitle(title.Data());
    hd->GetXaxis()->CenterTitle(1);
    hd->GetYaxis()->CenterTitle(1);

    // get range:
    double max = TMath::Max( hAtt[iatt]->hSummary->GetMaximum(),
                             hmc->GetMaximum());
    max = TMath::Max(  max, hd->GetMaximum());
    
    // draw histos
    hd->SetMaximum(max*1.2);
    hd->Draw();

    // don't draw buffer bins
    double xmax = hd->GetBinLowEdge(hd->GetNbinsX()) - (binbuffer*hd->GetBinWidth(1));
    double xmin = hd->GetBinLowEdge(1) - hd->GetBinWidth(1) + (binbuffer*hd->GetBinWidth(1));
    hd->GetXaxis()->SetRangeUser(xmin,xmax);
    hAtt[iatt]->hSummary->Draw("samee2");
    hmc->Draw("samehisto");
    hd->Draw("same");

  }

  //
  return;
}


/////////////////////////////////////////////////////////////
void fitPlots::drawFitSummary(int isamp, int ibin){

  // fill arrays
  fillArrays(ibin, isamp);

  // print at title
  TString title = Form("Sample %d Region %d",isamp,ibin+1);

  cc = new TCanvas("cc","cc",800,800);
  cc->Divide(2,2);
  TString xtitle[] = {"e/#mu PID",
                      "e/#pi^{0} PID",
                      "#mu/#pi PID",
                      "RC Parameter"};

  for (int iatt=0; iatt<4; iatt++){

    cc->cd(iatt+1);
    // calc summaryy
    hAtt[iatt]->calcSummary();
    hAtt[iatt]->hSummary->SetMarkerStyle(4);
    hAtt[iatt]->hSummary->GetXaxis()->SetTitle(xtitle[iatt].Data());
    hAtt[iatt]->hSummary->GetXaxis()->SetTitleSize(0.05);
    hAtt[iatt]->hSummary->GetXaxis()->SetNdivisions(5);
    hAtt[iatt]->hSummary->GetYaxis()->SetNdivisions(5);
    hAtt[iatt]->hSummary->SetTitle(title.Data());

    // draw MC original
    TH1D* hmc = hCompare->hManager->getSumHistogram(isamp,ibin,iatt,1);
    hmc->SetLineColor(kRed);

    // draw Data
    TH1D* hd = hCompare->hManager->hData[isamp][ibin][iatt];
    hd->SetMarkerStyle(8);

    // get range:
    double max = TMath::Max( hAtt[iatt]->hSummary->GetMaximum(),
                             hmc->GetMaximum());
    max = TMath::Max(  max, hd->GetMaximum());
    
    // draw
    hd->SetMaximum(max*1.2);
    hd->Draw();
    hAtt[iatt]->hSummary->Draw("samee2");
    hmc->Draw("samehisto");
    hd->Draw("same");
  }

  //
  return;
}



/////////////////////////////////////////////////////////////
void fitPlots::drawFitThrows(int isamp, int ibin){

  // fill arrays
  fillArrays(ibin, isamp);

  cc = new TCanvas("cc","cc",800,800);
  cc->Divide(2,2);

  for (int iatt=0; iatt<4; iatt++){

    cc->cd(iatt+1);

    // draw array
    hAtt[iatt]->drawArray();

    // draw MC original
    TH1D* hmc = hCompare->hManager->getSumHistogram(isamp,ibin,iatt,1);
    hmc->SetLineColor(kRed);
    hmc->Draw("samehisto");

    // draw Data
    TH1D* hd = hCompare->hManager->hData[isamp][ibin][iatt];
    hd->SetMarkerStyle(8);
    hd->Draw("samee");

  }


  return;
}

/////////////////////////////////////////////////////////////
void fitPlots::initArrays(){

  // seed arrays
  for (int iatt=0; iatt<nAtt; iatt++){
    TH1D* hseed = hCompare->hManager->getSumHistogramMod(0,5,iatt);
    hAtt[iatt] = new hArray(Form("att%d",iatt),hseed,nPoints);    
  }

  return;
}


/////////////////////////////////////////////////////////////
void fitPlots::fillArrays(int ibin, int isamp){
  
  // set rng
  TRandom2* rand = new TRandom2(nPoints);
  int nmax = parTree->GetEntries();

  // loop on points
  for (int ipt=0; ipt<nPoints; ipt++){
   
    // get random throw
    int ipar = rand->Integer(nmax);
    parTree->GetEntry(ipar);

    // apply pars
    applyPars();

    // fill arrays
    for (int iatt=0; iatt<nAtt; iatt++){
       
       // get mod histogram
       TH1D* hmod = hCompare->hManager->getSumHistogramMod(isamp, ibin, iatt);;
       
       // clone it in array
       hAtt[iatt]->setHistogram(ipt,hmod);

    } //< end att loop

  } //< end mcmc point loop

  return; 
}



/////////////////////////////////////////////////////////////
void fitPlots::applyPars(){

  // loop over mcmc pars and set int fitPars
  for (int ipar=0; ipar<mcmcPars->npars; ipar++){
   
    // get the index of the ipar-th parameter in the atmfitpars array
    int atmparindex = mcmcPars->parindex[ipar];

    // talk about it
    cout<<"set par "<<atmparindex<<" "<<" <- "<<ipar<<" "<<mcmcPars->par[ipar]<<endl;

    // change the atmfitpars array
    hCompare->thePars->setParameter(atmparindex, (double)mcmcPars->par[ipar]);
  }

  return;
}



#endif
