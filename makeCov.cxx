#include "TTree.h"
#include "TFile.h"
#include "TH2F.h"
#include "TMath.h"
#include "TLine.h"
#include "TString.h"
#include "TLatex.h"
#include "TCanvas.h"
#include <iostream>
#include "atmFitPars.h"

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

  //////////////////////
  //pull histograms
  TH1D* hpull;
  TH1D* hpulldist;

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
  void drawBinSubMatrix(int ibin);
//  void drawAlphaSubMatrix();
    
};


void makeCov::drawBinSubMatrix(int ibin){

}


void makeCov::drawLabeldCor(){
  TCanvas *cc = new TCanvas("cc","cc",700,700);
  cor->GetYaxis()->SetNdivisions(0);
  cor->GetXaxis()->SetNdivisions(0);
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

}
