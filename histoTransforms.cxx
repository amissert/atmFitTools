#ifndef HISTOTRANSFORMS_H
#define HISTOTRANSFORMS_H

#include "TH1D.h"
#include "TRandom3.h"
#include "TMath.h"
#include <time.h>
#include <iostream>
#include "TCanvas.h"

#ifndef GLOBAL_RANDOM
#define GLOBAL_RANDOM
TRandom3* randy = new TRandom3();
#endif

using namespace std;


/////////////////////////////////////////////////////
//returns a parameter for S/N ratio for a histogram
double getNoiseFactor(TH1D* hh){
//  if (hh->GetEntries()==0) return 0.;
  int nbins = hh->GetNbinsX();
  double S = 0.;
  double B = 0;
  for (int ibin = 1;ibin<=nbins;ibin++){
    double content = hh->GetBinContent(ibin);
    if (content>0.){
      S+=content;
      B++;
    }
  }
  if (B==0) return 0;
  S/=B;
  return 1./TMath::Sqrt(S);
}


///////////////////////////////////////////////////////////////////////////////////
//Custom smoothing method
void mySmooth(TH1D* hh,double factor=3.0){

  //////////////////////////////
  //get adjecent bin weights
  
  // use Gaussian weights
  double sigma = factor*getNoiseFactor(hh); //< set sigma using noise
  if (sigma==0) return; //no events in histogram
  double w2 = TMath::Gaus(2,0,sigma,1);
  double w1 = TMath::Gaus(1,0,sigma,1);
  double w0 = TMath::Gaus(0,0,sigma,1);

  ////////////////////////////////
  //clone in original histogram
  TH1D* htmp = (TH1D*)hh->Clone("tmphistosmooth");

  ///////////////////////////////////////////////
  //set new histogram contents from adjecent bins
  double newcontent;
  for (int ibin=0;ibin<hh->GetNbinsX();ibin++){
    newcontent = 0.;
    newcontent+=(htmp->GetBinContent(ibin-2)*w2);
    newcontent+=(htmp->GetBinContent(ibin-1)*w1);
    newcontent+=(htmp->GetBinContent(ibin)*w0);
    newcontent+=(htmp->GetBinContent(ibin+1)*w1);
    newcontent+=(htmp->GetBinContent(ibin+2)*w2);
    hh->SetBinContent(ibin,newcontent);
    hh->SetBinError(ibin,TMath::Sqrt(newcontent));
  }

  //normalize histogram
  double normscale = htmp->Integral()/hh->Integral();
  hh->Scale(normscale); 
  htmp->Delete(); 
  return; 
}

///////////////////////////////////
//Integral of box function
double B(double x,double a, double b){
  if (x<=a) return 0.;
  if (x>b) return 1.;
 // if (a==b) return 0;
  return (x-a)/(b-a);
}


//////////////////////////////////////////////////////////////////////////////////////
//useful test functions
TH1D* testBumpD(int nev,double sig=1.0,double mean=0.0,const char* name="testbumb"){
  TH1D* hset;
  hset->SetDefaultSumw2(kTRUE);
  TH1D* h = new TH1D(name,name,40,-10,10);
  for (int  i=0;i<nev;i++){
    h->Fill(randy->Gaus(mean,sig));
  }
  return h;
}

TH1D* testBump(int nev,double sig=1.0,double mean=0.0,const char* name="testbumb"){
  TH1D* hset;
  hset->SetDefaultSumw2(kTRUE);
  TH1D* h = new TH1D(name,name,40,-10,10);
  for (int  i=0;i<nev;i++){
    h->Fill(randy->Gaus(mean,sig));
  }
  return h;
}

TH1D* testTable(int nev,double width=4.0,double mean=0.0){
  TH1D* h = new TH1D("testbump","testbump",50,-10,10);
  for (int  i=0;i<nev;i++){
    h->Fill((randy->Rndm()*width) - (width/2.)+mean);
  }
  return h;
}
///////////////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////
//performs histogram modifications
void smearThisHisto(TH1D &hh, double spread, double bias=0.){

  //make sure the parameters are reasonable
  if (spread==0){
    cout<<"histoTransforms.cxx: smearThisHisto: cannont smear with 0 spread parameter!"<<endl;
    return;
  }

  //make temporary clone of input histogram 
  TH1D* htmp = (TH1D*)hh.Clone("htmp");

  //apply custom smooth function if statistics are low
//  if (hh.GetEntries()<10000.) mySmooth(htmp);

  //get some useful histogram parameters
  int nbins=hh.GetNbinsX();
  double binw = hh.GetBinWidth(1);

  //parameters for calculations
  double binedge;
  double sum;
  double binerr;
  double weight;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double smear = 1./spread;
  double mean = hh.GetMean() + (binw/2.);
  double shift = -1*(mean - (smear*mean)); //corrects for bias from smearing
  double sumw2;
  double sumw;


  //loop over bins and modify contents
  for (int newbin=1;newbin<=nbins;newbin++){

 //   if ((htmp->GetBinContent(newbin)<=0)||(htmp->GetBinContent(newbin-1)<=0)||(htmp->GetBinContent(newbin+1)<=0)){
 //     continue;
 //   }

    ////////////////////
    //set weight parameters
    sum = 0.;
    sumw = 0;
    sumw2 = 0;
    binerr=0.;
    binedge = htmp->GetBinLowEdge(newbin);
    ymin = ((binedge-bias)*smear) - shift;
    ymax = ymin + (binw*smear);
    for (int oldbin=1;oldbin<=nbins;oldbin++){
      xmin = htmp->GetBinLowEdge(oldbin);
      xmax = (xmin+binw);
      weight =  B(xmax,ymin,ymax)-B(xmin,ymin,ymax);
      double binc = htmp->GetBinContent(oldbin);
      sum+=(weight*htmp->GetBinContent(oldbin));
      binerr += weight*weight*(htmp->GetBinContent(oldbin)*htmp->GetBinContent(oldbin));
      sumw += weight;
    }
  //  if (sumw<=0.0001) continue;
    hh.SetBinContent(newbin,sum);
    hh.SetBinError(newbin,TMath::Sqrt(sum));
//    double scale = htmp->Integral()/hh.Integral();
//    hh.Scale(scale);
  }



  double scale = htmp->Integral()/hh.Integral();
  hh.Scale(scale);

  htmp->Delete();
  return;
}



//////////////////////////////////////////////////////////////
//Test function to compare to event-by-event
void compareEvtByEvt(int nevts, double mean0, double sig0,double scale, double bias){

  //nominal histogram
  TH1D* h0 = testBumpD(nevts, sig0, mean0, "nominal");
  double mean = h0->GetMean();

  //transformed histogram
  TH1D* htransform =  testBumpD(nevts, sig0, mean0, "transformed");
  htransform->SetLineColor(kRed);
  smearThisHisto(*htransform,scale,bias);

  //event by event histogram
  TH1D* hebe=  testBumpD(nevts, sig0, mean0, "evntbyevt");
  hebe->Reset();
  hebe->SetLineColor(kBlue);
  for (int i=0;i<nevts;i++){
    double xx = randy->Gaus(mean0,sig0); //< get random throw
    xx *= scale;
    xx += bias;
   // xx += bias*scale;
    xx -= scale*mean;
    xx += mean;
    hebe->Fill(xx);
  }

  TCanvas* cc = new TCanvas("cc","cc",800,700);
  h0->Draw();
  htransform->Draw("same");
  hebe->Draw("same");
  cc->Print("compare.png"); 

  return;

}



//smear it faster
void smearThisHistoFastMean(TH1D &hh, double* hcontent, double spread, double mean, double bias, double normscale=1.){

  //make sure the parameters are reasonable
  if (spread==0){
    cout<<"histoTransforms.cxx: smearThisHisto: cannont smear with 0 spread parameter!"<<endl;
    return;
  }

  //make temporary clone of input histogram 
//  TH1D* htmp = (TH1D*)hh.Clone("htmp");

  //apply custom smooth function if statistics are low
  // if (hh.GetEntries()<10000.) mySmooth(htmp);

  //get some useful histogram parameters
  int nbins=hh.GetNbinsX();
  double oldintegral = hh.Integral();
  double binw = hh.GetBinWidth(1);

  //parameters for calculations
  double binedge;
  double sum;
  double binerr;
  double weight;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double smear = 1./spread;
//  double mean = hh.GetMean() + (binw/2.);
  double shift = -1*(mean - (smear*mean)); //corrects for bias from smearing
  double sumw2;
  double sumw;

  //loop over bins and modify contents
  for (int newbin=1;newbin<=nbins;newbin++){
    sum = 0.;
    sumw = 0;
    sumw2 = 0;
    binerr=0.;
    binedge = hh.GetBinLowEdge(newbin);
    ymin = ((binedge-bias)*smear) - shift;
    ymax = ymin + (binw*smear);
    for (int oldbin=1;oldbin<=nbins;oldbin++){
      xmin = hh.GetBinLowEdge(oldbin);
      xmax = (xmin+binw);
      weight =  B(xmax,ymin,ymax)-B(xmin,ymin,ymax);
      double binc = hcontent[oldbin];
      sum+=(weight*binc);
      binerr += weight*weight*(binc*binc);
      sumw += weight;
    }
    hh.SetBinContent(newbin,sum);
  }
  double newintegral = hh.Integral();
  double scale = oldintegral/newintegral;
  hh.Scale(scale*normscale);
  return;
}



//smear it faster
void smearThisHistoFast(TH1D &hh, double* hcontent, double spread, double bias=0.){

  if (hh.Integral()<1e-4) return;
  //make sure the parameters are reasonable
  if (spread==0){
    cout<<"histoTransforms.cxx: smearThisHisto: cannont smear with 0 spread parameter!"<<endl;
    return;
  }

  //get some useful histogram parameters
  int nbins=hh.GetNbinsX();
  double oldintegral = hh.Integral();
  double binw = hh.GetBinWidth(1);

  //parameters for calculations
  double binedge;
  double sum;
  double binerr;
  double weight;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double smear = 1./spread;
  double mean = hh.GetMean() + (binw/2.);

  double shift = -1*(mean - (smear*mean)); //corrects for bias from smearing
  double sumw2;
  double sumw;

  //loop over bins and modify contents
  for (int newbin=1;newbin<=nbins;newbin++){
    sum = 0.;
    sumw = 0;
    sumw2 = 0;
    binerr=0.;
    binedge = hh.GetBinLowEdge(newbin);
    ymin = ((binedge-bias)*smear) - shift;
    ymax = ymin + (binw*smear);
    for (int oldbin=1;oldbin<=nbins;oldbin++){
      xmin = hh.GetBinLowEdge(oldbin);
      xmax = (xmin+binw);
      weight =  B(xmax,ymin,ymax)-B(xmin,ymin,ymax);
      double binc = hcontent[oldbin];
      sum+=(weight*binc);
      binerr += weight*weight*(binc*binc);
      sumw += weight;
    }
    hh.SetBinContent(newbin,sum/sumw);
  }
  double newintegral = hh.Integral();
  double scale = oldintegral/newintegral;
  hh.Scale(scale);
  return;
}


double testtime(){
  int ntry = 25000;
  clock_t t1,t2;
  TH1D* hb = testBump(1000);
  TH1D* hmod = (TH1D*)hb->Clone("htmp");
  t1=clock();
  for (int i=0;i<ntry;i++){
    smearThisHisto(*hb,1.1,1.2);
  }
  t2=clock();
  double diff = ((double)t2-(double)t1)/((double)ntry);
  cout<<"time: "<<diff<<endl;
  return diff;
  
}




#endif

