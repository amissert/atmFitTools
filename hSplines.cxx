#ifndef HSPLINES_C
#define HSPLINES_C

#include "hSplines.h"

using namespace std;

void hSplines::initFlgs(){
  for (int i=0; i<nHistoBins; i++){
    for (int j=0; j<nSyst; j++){
      setflg[i][j] = 0;
    }
  }
  return;
}

double hSplines::evaluateSpline(int ibin, int ipar, double parvalue){
#ifndef T2K
  if (parvalue<0.) parvalue = 0.;
#endif
  if (theSpline[ibin][ipar]==NULL) return 1.0;
  return (double)theSpline[ibin][ipar]->Eval(parvalue);
}

void  hSplines::drawSpline(int ibin, int isyst){
  if (theSpline[ibin][isyst]==NULL){
    cout<<"no such spine exists!"<<endl;
    return;
  }
  theSpline[ibin][isyst]->Draw();
}

void  hSplines::draw2D(int npts,int isyst){
  cout<<"getting par max and min"<<endl;
  double parmin= theSpline[1][isyst]->GetXmin();
  double parmax= theSpline[1][isyst]->GetXmax(); 
  double xx;
  double parval;
  double increment = (parmax - parmin) / npts;
  cout<<"delete prev histo"<<endl;
//  if (drawHisto!=NULL) drawHisto->Delete();
  cout<<"make new histo histo"<<endl;
  drawHisto= new TH2F("hdraw","hdraw",nHistoBins,baseHisto->GetBinLowEdge(1),
                      (baseHisto->GetBinLowEdge(baseHisto->GetNbinsX()) + baseHisto->GetBinWidth(baseHisto->GetNbinsX())),
                      npts,parmin-increment/2.,parmax+increment/2.);
  cout<<"fill xy histogram"<<endl;
  cout<<"nbins"<<nHistoBins<<endl;
  for (int xbin=1;xbin<=drawHisto->GetNbinsX();xbin++){
    for (int ybin=1;ybin<drawHisto->GetNbinsY();ybin++){
      xx = baseHisto->GetXaxis()->GetBinCenter(xbin); 
      parval = drawHisto->GetYaxis()->GetBinCenter(ybin);
      cout<<"parameter value: "<<parval<<endl;
      cout<<"spline value: "<< theSpline[xbin-1][isyst]->Eval(parval)<<endl;
      // old way 
//      drawHisto->SetBinContent(xbin,ybin,theSpline[xbin-1][isyst]->Eval(parval)*baseHisto->GetBinContent(xbin)); 
      // new way
      double binc = evaluateSpline(xbin-1,isyst,parval)*baseHisto->GetBinContent(xbin);
      drawHisto->SetBinContent(xbin,ybin,binc); 
    }
  }
  cout<<"draw histogram"<<endl;
  drawHisto->SetContour(20);
  drawHisto->Draw("lego2");
  return;
}

void  hSplines::debugTest(){
  //test it
  double x1[5] = {0.,1.,2.,3.,4.};
  double y1[5] = {0.5,0.1,0.9,1.3,1.1};
  nSyst = 1;
  nHistoBins = 5;
  baseHisto = new TH1D("hdebug","hdebug",nHistoBins,0,5);
  baseHisto->Fill(1.);
  baseHisto->Fill(1.);
  baseHisto->Fill(2.);
  baseHisto->Fill(3.);
  baseHisto->Fill(4.);
  baseHisto->Fill(5.);
  baseHisto->Fill(0.);
  baseHisto->Fill(4.);
  baseHisto->Fill(4.);
  baseHisto->Fill(4.);
  for (int ibin=1;ibin<=nHistoBins;ibin++){
    buildSpline(ibin,0,x1,y1,5);
  }
  draw2D(10,0);
  return;
}

void hSplines::setSpline(int ibin, int isyst, TSpline3 *spline){
  theSpline[ibin][isyst] = spline;
  setflg[ibin][isyst] = 1;
  checkSum--;
  return;
}

void hSplines::buildSpline(int ibin, int isyst,double* X, double*Y, int N){
  //temporary arrays
  const int npts = N;
  double    xvals[npts];
  double    yvals[npts];
  //ignore x values less than zero:
  int index=0;
  for (int ipt=0;ipt<npts;ipt++){
#ifdef T2K
    if (Y[ipt]<0) std::cout<<"!!!!! WEIGHT < 0 !!!!!"<<std::endl;
    xvals[index]=X[ipt];
    yvals[index]=Y[ipt];
    index++;
#else
    if (Y[ipt]>=0.){
      xvals[index]=X[ipt];
      yvals[index]=Y[ipt];
      index++;
    }
#endif
  }
  TString splineName = nameTag.Data();
  splineName.Append(Form("_spline_bin%d_par%d",ibin,isyst));
 // if (!theSpline[ibin][isyst])  theSpline[ibin][isyst] = new TSpline3(splineName.Data(),xvals,yvals,index); 
  theSpline[ibin][isyst] = new TSpline3(splineName.Data(),xvals,yvals,index); 
  theSpline[ibin][isyst]->SetName(splineName.Data());
  cout<<"hSplines: Built spline: "<<splineName.Data()<<endl;
  setflg[ibin][isyst] = 1;  
  checkSum--;
  return;
}

hSplines::hSplines(TH1D* h, int nsyst, const char* name){

   // fix the name of this histogram + spline object
   if (name!=NULL) nameTag = name;
   else{
     nameTag = h->GetName();
   }

   // get number of histo bins and number of systematic parameters
   nHistoBins=h->GetNbinsX();
   nSyst=nsyst;

   //initialize base histogram
   baseHisto=h;

   //initialize modified histogram
   TString modhname = baseHisto->GetName();
   modhname.Append("_modified");
//   int hbins = baseHisto->GetNbinsX();
 //  double hmin = baseHisto->GetBinLowEdge(1);
  // double hmax = baseHisto->GetBinLowEdge(baseHisto->GetNbinsX())+baseHisto->GetBinWidth(1);
   //modHisto = new TH1D(modhname.Data(),modhname.Data(),hbins,hmin,hmax);
   modHisto = (TH1D*)baseHisto->Clone(modhname.Data());
   modHisto->Reset();
   cout<<"mod histogram: "<<modHisto->GetName()<<endl;

 //  int basenbins = baseHisto->GetNbinsX();
 //  double basexmin  = baseHisto->GetBinLowEdge(1);
  // double  basexmax = baseHisto->GetBinLowEdge(basenbins);
   checkSum = nHistoBins*nSyst;
   initFlgs(); //< initialize set spline flag values
   //modHisto= new TH1D(modhname.Data(),modhname.Data(),basenbins,basexmin,basexmax);
   
//   modHisto=(TH1D*)baseHisto->Clone(modhname.Data());
//   modHisto->Reset();
}

#endif
