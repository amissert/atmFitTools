#ifndef FITPLOTS_C
#define FITPLOTS_C

#include "fitPlots.h"
#include "hBound.cxx"



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
void fitPlots::drawFitSummary(int isamp, int ibin){

  // fill arrays
  fillArrays(ibin, isamp);

  // print at title
  TString title = Form("Sample %d Region %d",isamp,ibin+1);
//  TLatex* tit = new TLatex(0.5,0.95,title.Data());
//  tit->SetTextAlign(22);

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
