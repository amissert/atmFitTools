#include "TTree.h"
#include "TFile.h"
#include "TH2F.h"
#include "TMath.h"
#include "TLine.h"
#include "TF1.h"
#include "TString.h"
#include "TLatex.h"
#include "TCanvas.h"
#include <iostream>
#include <vector>
#include "atmFitPars.h"
#include "parMap.C"

#define NTOTALPARS 500

using namespace std;


/////////////////////////////////////////////////////////
//make covariance on correlation matrix from markov chain
class makeCov{

  public:

  ////////////////
  //constructor
  makeCov(const char* parfile = "");

  /////////////////
  //parameter file
  TString runparfile;

  /////////////////
  //tree of mcmc steps
  TTree* partree;
  // tree leaves
  float par[500];
  float parnominal[500];
  int parbin[500];
  int parcomp[500];
  int paratt[500];
  int parsyst[500];
  int parindex[500];
  int npar;


  ////////////////////////////////////////
  // parameter labels
  parMap* pmap;

  ////////////////////////////////////////
  //set parametet tree
  void setParTree(TTree* tr);

  ///////////////////////////////////////////
  //set the number of burn in steps to ignore
  int nburn;

  ///////////////////////////////////////////
  //set the number of steps to skip 
  int nskip;

  /////////////////////
  //matrix histograms
  TH2F* cov;
  TH2F* cor;
  TH2D* hsubcor;

  //////////////////////
  //pull histograms
  TH1D* hpull;
  TH1D* hpulldist;
  TH1D* hvalues;
  TH1D* hsysprior;
  TH1D* hsysvalues;

  ///////////////////////
  //dividing lines
  TLine* lhoriz;
  TLine* lvert;

  ////////////////////////
  // number of parameters
  int nsyspar;
  int nfitpartot;

  ////////////////////////
  //arrays
  double parmean[NTOTALPARS];
  double parsigma[NTOTALPARS];
  double pardefault[NTOTALPARS];
  double covarray[NTOTALPARS][NTOTALPARS];
  double corarray[NTOTALPARS][NTOTALPARS];

  /////////////////////////
  //build that matrix
  void buildMatrix();
  TH2D* getBinSubMatrix(int ibin, int itype=1);
  TH2D* getAlphaSubMatrix(int nsub=2);

  //////////////////////////////
  // draw nice correlation matrix
  void drawCor();

  ////////////////////////////////
  // print all 1d parameter distributions to a directory
  void printall1D(const char* dir);

  ////////////////////////////////////
  // print some par info
  void printParInfo(int ipar);

  ///////////////////////////////////
  // print errors
  void printParErrors();

  ///////////////////////////////////
  // draw with lines and labels
  void drawLabeldCor();
  void drawVertBinLines();
  void drawHorzBinLines();
  void drawBinSubMatrix(int ibin, int itype=1);
  void drawValues(int icomp, int iatt, int itype);
  void drawSysValues();

  void setSysPrior();
  void drawAlphaSubMatrix(int nbin);
    
  TF1* zero;
  TF1* fone;

};


///////////////////////////////////////////////////// 
void makeCov::setSysPrior(){
  float labelsize = 0.05;
  hsysprior = new TH1D("hsysprior","hsysprior",37,0,37);
  hsysprior->SetFillColor(kGray);
  for (int ibin=1; ibin<=37; ibin++){
    hsysprior->SetBinContent(ibin,1);
    hsysprior->GetXaxis()->SetBinLabel(ibin,pmap->getSystParName(ibin-1));
  }
  hsysprior->GetXaxis()->LabelsOption("v");
  hsysprior->GetXaxis()->SetLabelSize(labelsize);
  hsysprior->SetBinError(1,1.0);
  hsysprior->SetBinError(2,.411);
  hsysprior->SetBinError(3,.216);
  hsysprior->SetBinError(4,.155);
  hsysprior->SetBinError(5,.125);
  hsysprior->SetBinError(6,.105);
  hsysprior->SetBinError(7,.0805);
  hsysprior->SetBinError(8,.066);
  hsysprior->SetBinError(9,.0542);
  hsysprior->SetBinError(10,.0398);
  hsysprior->SetBinError(11,.0344);
  hsysprior->SetBinError(12,.0226);
  hsysprior->SetBinError(13,.0165);
  hsysprior->SetBinError(14,.009);
  hsysprior->SetBinError(15,.25);
  hsysprior->SetBinError(16,.15);
  hsysprior->SetBinError(17,.2);
  hsysprior->SetBinError(18,.2);
  hsysprior->SetBinError(19,.05);
  hsysprior->SetBinError(20,.1);
  hsysprior->SetBinError(21,.1);
  hsysprior->SetBinError(22,.1);
  hsysprior->SetBinError(23,.1);
  hsysprior->SetBinError(24,.1);
  hsysprior->SetBinError(25,.1);
  hsysprior->SetBinError(26,.1);
  hsysprior->SetBinError(27,.1);
  hsysprior->SetBinError(28,.1);
  hsysprior->SetBinError(29,.1);
  hsysprior->SetBinError(30,.1);
  hsysprior->SetBinError(31,.1);
  hsysprior->SetBinError(32,.1);
  hsysprior->SetBinError(33,.1);
  hsysprior->SetBinError(34,.1);
  hsysprior->SetBinError(35,.1);
  hsysprior->SetBinError(36,.1);
  hsysprior->SetBinError(37,.1);
  hsysprior->SetBinError(38,.1);

  return;
}



///////////////////////////////////////////////////// 
void makeCov::drawSysValues(){

  TCanvas* cc = new TCanvas("cc","cc",1400,700);
  cc->GetPad(0)->SetBottomMargin(0.4);

  int nsyspar = 37;
  if (hsysvalues!=NULL){
    hsysvalues->Delete();
  }
  hsysvalues = new TH1D("hsysval","hsysval",nsyspar,0,nsyspar);
  hsysprior->Draw("e2");
  for (int isys=0; isys<nsyspar; isys++){
    int index = pmap->getSysIndex(isys);
    hsysvalues->SetBinContent(isys+1,parmean[index]);    
    hsysvalues->SetBinError(isys+1,parsigma[index]);    
    hsysvalues->GetXaxis()->SetBinLabel(isys+1,pmap->getSystParName(isys).Data());
  }
 
  hsysvalues->GetXaxis()->LabelsOption("v");
  hsysvalues->GetYaxis()->SetTitle("value");
  hsysvalues->SetLineColor(kRed);
  hsysvalues->SetLineWidth(3);
  hsysvalues->SetMarkerStyle(1);
  hsysvalues->Draw("same");
  return;
}


///////////////////////////////////////////////////// 
void makeCov::drawValues(int icomp, int iatt, int itype){

  TCanvas* cc = new TCanvas("cc","cc",700,900);
  cc->GetPad(0)->SetBottomMargin(0.3);
  int nbins = 6;

  if (hvalues!=NULL){
    hvalues->Delete();
  }

  hvalues = new TH1D("hvalues","hvalues",nbins,0,nbins);

  // fill values
  for (int i=0; i<nbins; i++){
     int index = pmap->getGlobalIndex(i,icomp,iatt,itype);
     hvalues->SetBinContent(i+1,parmean[index]);
     hvalues->SetBinError(i+1,parsigma[index]);
     hvalues->GetXaxis()->SetBinLabel(i+1,pmap->getName(index).Data());
  }

  double absmax = TMath::Max(TMath::Abs(hvalues->GetMinimum()),TMath::Abs(hvalues->GetMaximum()));
  hvalues->SetMaximum(2.*(absmax - (float)(itype+1)));
  hvalues->SetMinimum(-2.*(absmax - (float)(itype+1)));
  hvalues->GetXaxis()->LabelsOption("v");
  hvalues->GetYaxis()->SetTitle("Value");
  hvalues->SetLineColor(kRed);
  hvalues->SetLineWidth(3);
  hvalues->SetLabelSize(0.04);
  hvalues->Draw();
 
  return;
}


///////////////////////////////////////////////////// 
void makeCov::drawAlphaSubMatrix(int nsub){

  TCanvas *cc = new TCanvas("cc","cc",700,750);
  cc->GetPad(0)->SetLeftMargin(0.25);
  cc->GetPad(0)->SetRightMargin(0.12);
  cc->GetPad(0)->SetBottomMargin(0.25);
  getAlphaSubMatrix(nsub);
  hsubcor->Draw("colz");
  return;
}


///////////////////////////////////////////////////// 
void makeCov::drawBinSubMatrix(int ibin, int itype){

  TCanvas *cc = new TCanvas("cc","cc",800,700);
  cc->GetPad(0)->SetLeftMargin(0.35);
  cc->GetPad(0)->SetRightMargin(0.12);
  cc->GetPad(0)->SetBottomMargin(0.35);
  getBinSubMatrix(ibin,itype);

  // lines to separate blocks
  double lmin  = 0;
  double lmax  = (double)hsubcor->GetNbinsX();
  TLine* lv[3];
  TLine* lh[3];
  double sep = 6;
  for (int i=0; i<3; i++){
    lv[i] = new TLine(sep,lmin,sep,lmax);
    lh[i] = new TLine(lmin,sep,lmax,sep);
    sep += 6;
  }
  hsubcor->SetContour(100);
  hsubcor->SetBit(14,1);
  hsubcor->Draw("colz");
  for (int i=0; i<3; i++){
    lv[i]->Draw("same");
    lh[i]->Draw("same");
  }
  return;
}


///////////////////////////////////////////////////// 
TH2D* makeCov::getAlphaSubMatrix(int nsub){

  // number of things
  int nsubins = 0;
  if (nsub==0) nsubins = 19; //< flux and xsec
  if (nsub==1) nsubins = 18; //< normalization
  if (nsub==2) nsubins = 37; //< both

  // delete if exists
  if (hsubcor!=NULL){
    hsubcor->Delete();
  }

  // make new
  hsubcor = new TH2D("hsubcor","hsubcor",nsubins,0,nsubins,nsubins,0,nsubins);

  // get indicies
  int istart = 0;
  if (nsub==0) istart = pmap->nParsTot-37;
  if (nsub==1) istart = pmap->nParsTot-18;
  if (nsub==2) istart = pmap->nParsTot-37;

  // get values
  for (int ii=0; ii<nsubins; ii++){
    for (int jj=0; jj<nsubins; jj++){
//      if (ii==jj) continue;
      hsubcor->SetBinContent(ii+1,jj+1,corarray[istart+ii][istart+jj]);
    }
    hsubcor->GetXaxis()->SetBinLabel(ii+1,pmap->getName(istart+ii));
    hsubcor->GetYaxis()->SetBinLabel(ii+1,pmap->getName(istart+ii));
  }

  // set labels
  hsubcor->GetXaxis()->LabelsOption("v");
  double maxcor = TMath::Max(TMath::Abs(hsubcor->GetMinimum()),
                             TMath::Abs(hsubcor->GetMaximum()));
  hsubcor->SetMinimum(-1*maxcor);
  hsubcor->SetMaximum(maxcor);
  hsubcor->SetContour(100);

  return hsubcor;
}




///////////////////////////////////////////////////// 
TH2D* makeCov::getBinSubMatrix(int ibin, int itype){

  // number of things
  int ncomp = 6;
  int natt = 4;
  int ntype = 1;
  int nsubins = ncomp*natt*ntype;
  double labsize = 0.035;

  // delete if exists
  if (hsubcor!=NULL){
    hsubcor->Delete();
  }

  // make new
  hsubcor = new TH2D("hsubcor","hsubcor",nsubins,0,nsubins,nsubins,0,nsubins);

  // get indicies
  int istart = pmap->getGlobalIndex(ibin,0,0,0);
  int iend = istart + nsubins;

  // get bins
  vector<int> vbins;
    for (int iatt=0; iatt<natt; iatt++){
      for (int icomp=0; icomp<ncomp; icomp++){
//  for (int icomp=0; icomp<ncomp; icomp++){
//    for (int iatt=0; iatt<natt; iatt++){
      int binum = pmap->getGlobalIndex(ibin,icomp,iatt,itype);
      vbins.push_back(binum);
    }
  }

  // get values
  for (int ii=0; ii<nsubins; ii++){
    for (int jj=0; jj<nsubins; jj++){
//      if (ii==jj) continue;
      hsubcor->SetBinContent(ii+1,jj+1,corarray[vbins.at(ii)][vbins.at(jj)]);
    }
    hsubcor->GetXaxis()->SetBinLabel(ii+1,pmap->getName(vbins.at(ii)).Data());
    hsubcor->GetYaxis()->SetBinLabel(ii+1,pmap->getName(vbins.at(ii)).Data());
  }

  // set labels
  hsubcor->GetXaxis()->LabelsOption("v");

  double maxcor = TMath::Max(TMath::Abs(hsubcor->GetMinimum()),
                             TMath::Abs(hsubcor->GetMaximum()));
  hsubcor->SetMinimum(-1*maxcor);
  hsubcor->SetMaximum(maxcor);
  hsubcor->GetXaxis()->SetLabelSize(labsize);
  hsubcor->GetYaxis()->SetLabelSize(labsize);

  return hsubcor;
}



void makeCov::drawLabeldCor(){
  TCanvas *cc = new TCanvas("cc","cc",700,700);
  cor->GetYaxis()->SetNdivisions(0);
  cor->GetXaxis()->SetNdivisions(0);
  cor->SetMinimum(-1.);
  cor->SetMaximum(1.);
  cor->Draw("colz");
  drawVertBinLines();
  drawHorzBinLines();
  return;
}



void makeCov::drawHorzBinLines(){

 TLine* L[6];
 int spacing = 48;
 int offset  = -20;
 int max     = 325;
 for (int i=0; i<6; i++){
    double y = spacing + (double)i*spacing;
    double center = y - (spacing)/2.;
    TString label = Form("Region %d",i+1);
    TLatex* texlabel = new TLatex(offset,center,label.Data());
    texlabel->SetTextAlign(22);
    texlabel->SetTextAngle(90);
    texlabel->SetTextSize(0.025);
    texlabel->Draw();
    L[i] = new TLine(offset,y,max,y);
    L[i]->Draw();
 }

 double alphacenter = (325.-288)/2. + 288.;
 TString label = "#alpha";
 TLatex* texlabel = new TLatex(offset,alphacenter,label.Data());
 texlabel->SetTextAlign(22);
 texlabel->SetTextAngle(90);
 texlabel->SetTextSize(0.05);
 texlabel->Draw();

 return;
}



void makeCov::drawVertBinLines(){

 TLine* L[6];
 double spacing = 48;
 double offset  = -20;
 double max     = 325;
 for (int i=0; i<6; i++){
    double x = spacing + (double)i*spacing;
    double center = x - (spacing)/2.;
    TString label = Form("Region %d",i+1);
    TLatex* texlabel = new TLatex(center,offset,label.Data());
    texlabel->SetTextAlign(22);
    texlabel->SetTextSize(0.025);
    texlabel->Draw();
    L[i] = new TLine(x,offset,x,max);
    L[i]->Draw();
 }

 double alphacenter = (325.-288)/2. + 288.;
 TString label = "#alpha";
 TLatex* texlabel = new TLatex(alphacenter,offset,label.Data());
 texlabel->SetTextAlign(22);
// texlabel->SetTextAngle(90);
 texlabel->SetTextSize(0.05);
 texlabel->Draw();
 return;
}


void makeCov::printParErrors(){
  for (int ipar=0; ipar<nfitpartot; ipar++){
    double fiterr = parsigma[ipar];
    double shifterr = TMath::Abs(parmean[ipar] - pardefault[ipar]);
    cout<<"Par "<<ipar<<" Error: "<<fiterr+shifterr;
    cout<<" ("<<shifterr<<" shift and "<<fiterr<<" fit)"<<endl;
  }
  return;
}


void makeCov::printParInfo(int ipar){

     cout<<"par bin: "<<parbin[ipar]<<endl;
     cout<<"par comp: "<<parcomp[ipar]<<endl;
     cout<<"par att: "<<paratt[ipar]<<endl;

     return;
}


void makeCov::setParTree(TTree* tr){
  //setup mcmc trees
  partree = tr;
  partree->SetBranchAddress("par",par);
  partree->SetBranchAddress("npars",&npar);
  partree->SetBranchAddress("parnominal",parnominal);
  partree->SetBranchAddress("parbin",parbin);
  partree->SetBranchAddress("paratt",paratt);
  partree->SetBranchAddress("parcomp",parcomp);
  partree->SetBranchAddress("parindex",parindex);

  partree->GetEntry(0); //fills npar
  cout<<"Total # of parameters: "<<npar<<endl;
  cout<<"Total # of steps: "<<partree->GetEntries()<<endl;
  cout<<"Burn-in: "<<nburn<<endl;

  return;

}

void makeCov::printall1D(const char* dir){

  TString dirname = dir;

  TString basename = "h1D_parameter";
  TString branchname = "par";
;
  int npts = partree->GetEntries();

  //setup mcmc trees
/*  double par[500];
  double parnominal[500];
  int npar;
  partree->SetBranchAddress("par",par);
  partree->SetBranchAddress("npars",&npar);
  partree->SetBranchAddress("parnominal",parnominal);
  partree->GetEntry(0); //fills npar
  cout<<"Total # of parameters: "<<npar<<endl;
  cout<<"Total # of steps: "<<partree->GetEntries()<<endl;
  cout<<"Burn-in: "<<nburn<<endl;
*/
  // canvas setup
  TCanvas* cc = new TCanvas("cc","cc",800,700);
  //partree->Draw("par[0]");
 // cc->Print("testplot.png");

  // print
  for (int ipar=0; ipar<npar; ipar++){
    TString bname = branchname.Data();
    bname.Append(Form("[%d]",ipar));
    partree->Draw(bname.Data());
    TString plotname = dirname.Data();
    plotname.Append(basename.Data());
    plotname.Append(Form("%d.png",ipar));
    cc->Print(plotname.Data());
  }

  return;      
}

void makeCov::drawCor(){
  
  cor->SetStats(0);
  cor->GetZaxis()->SetRangeUser(-1.,1.);
  cor->Draw("colz");

  double xmax = cor->GetXaxis()->GetXmax();
  double ymax = cor->GetYaxis()->GetXmax();
  double xmin = cor->GetXaxis()->GetXmin();
  double ymin = cor->GetYaxis()->GetXmin();
  double xsep = (double)nfitpartot - (double)nsyspar;
  lhoriz = new TLine(0,xsep,xmax,xsep);
  lvert = new TLine(xsep,0,xsep,ymax);
  lvert->SetLineWidth(3);
  lhoriz->SetLineWidth(3);
  lhoriz->Draw("same");
  lvert->Draw("same");

  return;
}



void makeCov::buildMatrix(){

  //create matrix templates
  cov = new TH2F("cov","cov",npar,0.,(double)npar,npar,0.,(double)npar);
  cor = new TH2F("cor","cor",npar,0.,(double)npar,npar,0.,(double)npar);

  const int npartot = npar; //< total number of fit parameters in this mcmc cloud
  nfitpartot = npar;

  //calc means
  int npts = partree->GetEntries()-nburn;
  double sumpts=0;
  int nskipped=0;

  // loop over mcmc opints
  for (int iev=nburn; iev<partree->GetEntries(); iev++){
    nskipped++;
    if (nskipped<nskip) continue;
    nskipped = 0;
    partree->GetEntry(iev);  //< read in parameters from MCMC cloud
    // loop over fit parameters
    for (int ipar=0; ipar<npartot; ipar++){
      parmean[ipar]+=(par[ipar]);
    }
    sumpts++;
  }

  // normalize
  for (int ipar=0; ipar<npartot; ipar++){
    parmean[ipar]/=sumpts;
  }

  // print
  for (int kk=0;kk<npartot;kk++){
      cout<<"mean: "<<kk<<" "<<parmean[kk]<<endl;
  }

  //calc matrix
  double norm = 1./(sumpts-1.);
  nskipped=0;
  for (int jev=nburn; jev<partree->GetEntries(); jev++){
    nskipped++;
    if (nskipped<nskip) continue;
    nskipped = 0;
    partree->GetEntry(jev);
    for (int i0=0;i0<npartot;i0++){
      for (int j0=0;j0<npartot;j0++){
        covarray[i0][j0] += ((norm)*( (par[i0]-parmean[i0]) * (par[j0]-parmean[j0]) ));
      }
    }
  }
  for (int kk=0;kk<npartot;kk++){
      cout<<"matrix: "<<kk<<" "<<covarray[kk][kk]<<endl;
  };

  //fill histogram of matrix values
  for (int j=0;j<npartot;j++){
    for (int k=0;k<npartot;k++){
      cov->SetBinContent(j+1,k+1,covarray[j][k]);
      cor->SetBinContent(j+1,k+1,((covarray[j][k])/sqrt((covarray[j][j]*covarray[k][k]))));
      corarray[j][k] =  ((covarray[j][k])/sqrt( (covarray[j][j]*covarray[k][k]) ));
    }
    parsigma[j] = TMath::Sqrt( covarray[j][j] );
  }

  nsyspar=0;
  partree->GetEntry(1);
  for (int ipar=0; ipar<npartot; ipar++){
    pardefault[ipar] = parnominal[ipar];
  }

  // make pull histogram
  hpull = new TH1D("hpull","hpull",npartot,0,npartot);
  hpulldist = new TH1D("hpulldist","hpulldist",25,-5,5);
  for (int ipar=0; ipar<npartot; ipar++){
    double pullvalue = (parmean[ipar] - pardefault[ipar])/parsigma[ipar];
    hpull->Fill(ipar, pullvalue);
    hpulldist->Fill(pullvalue);
  }

  // remove pull errors
  for (int ibin=1; ibin<hpull->GetNbinsX(); ibin++){
    hpull->SetBinError(ibin,0.);
  }

  // draw the correlation
  cor->SetContour(100);
  cor->Draw("colz");
  
  ///////////////////////
  return;

}

makeCov::makeCov(const char* parfile){

  runparfile = parfile;

 // fix initial array values to zero
 for (int i=0; i<NTOTALPARS; i++){
   parmean[i]=0.;
   parsigma[i]=0.;
   pardefault[i]=0.;
   for (int j=0; j<NTOTALPARS; j++){
      corarray[i][j] = 0.;
      covarray[i][j] = 0.;
   }
 }

 pmap = new parMap();
 zero = new TF1("zer","0",-1000,1000);
 fone = new TF1("one","1",-1000,1000);
 setSysPrior();
}



