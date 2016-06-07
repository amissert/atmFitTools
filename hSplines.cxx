
  
#include "hSplines.h"

double hSplines::evaluateSpline(int ibin, int ipar, double parvalue){
  //if (parvalue<0.) parvalue = 0.;
  if (theSpline[ibin][ipar]==NULL) return 1.0;
  else return (double)theSpline[ibin][ipar]->Eval(parvalue);
}

void  hSplines::drawSpline(int ibin, int isyst){
  if (theSpline[ibin][isyst]==NULL){
    cout<<"no such spine exists!"<<endl;
    return;
  }
  theSpline[ibin][isyst]->Draw();
}

void  hSplines::draw2D(int npts,int isyst){
  double parmin= theSpline[1][isyst]->GetXmin();
  double parmax= theSpline[1][isyst]->GetXmax(); 
  cout<<npts<<", getting par max and min "<<parmax<<" "<<parmin<<endl;
  double increment = (parmax-parmin)/npts;
  double xx;
  double parval;
  //  cout<<"delete prev histo"<<endl;
  if (drawHisto!=NULL) drawHisto->Delete();
  //cout<<"make new histo histo"<<endl;
  drawHisto= new TH2F("hdraw","hdraw",nHistoBins,baseHisto->GetBinLowEdge(1),
                      (baseHisto->GetBinLowEdge(baseHisto->GetNbinsX()) + baseHisto->GetBinWidth(baseHisto->GetNbinsX())),
                      npts,parmin-increment/2.,parmax+increment/2.);
  cout<<"fill xy histogram"<<endl;
  cout<<"nbins"<<nHistoBins<<endl;
  for (int xbin=1;xbin<=drawHisto->GetNbinsX();xbin++){
    for (int ybin=1;ybin<=drawHisto->GetNbinsY();ybin++){
      xx = baseHisto->GetXaxis()->GetBinCenter(xbin); 
      parval = drawHisto->GetYaxis()->GetBinCenter(ybin);
      if (xbin==1) {
	cout<<xbin<<", parameter value: "<<parval<<endl;
	//cout<<"spline value: "<< theSpline[xbin][isyst]->Eval(parval)<<endl;
      }
      drawHisto->SetBinContent(xbin,ybin,theSpline[xbin][isyst]->Eval(parval)); 
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

void hSplines::buildModHisto(int ipar, double parval){ 

  
  double newcontent;
  for (int ibin=1;ibin<=nHistoBins;ibin++){
    newcontent = baseHisto->GetBinContent(ibin);
    //if (theSpline[ibin][ipar]->Eval(parval)==0){
    //   cout<<"zero spline interpolation for bin"<<ibin<<" with content "<<newcontent<<endl;
    //}
    newcontent *= theSpline[ibin][ipar]->Eval(parval); 
    modHisto->SetBinContent(ibin,newcontent);
  }
  return;
}

TH1D* hSplines::buildModHistoAllPar(int npars, double *systPars){
  double binx;
  double newcontent;
  double oldcontent;
  double weightsum;

    TString modhname = baseHisto->GetName();
    modhname.Append("_modified");
    modHisto=(TH1D*)baseHisto->Clone(modhname.Data());
  
  for (int ibin=1;ibin<=nHistoBins;ibin++){
    newcontent = baseHisto->GetBinContent(ibin);
    binx = baseHisto->GetBinCenter(ibin);
    weightsum=1.;
    for (int isyst=0;isyst<npars;isyst++){
      weightsum *= theSpline[ibin][isyst]->Eval(systPars[isyst]);
    }
    newcontent*=weightsum;
    modHisto->SetBinContent(ibin,newcontent);
  }
  return modHisto;
}

void hSplines::setSpline(int ibin, int isyst, TSpline3 *spline){
  theSpline[ibin][isyst] = spline;
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
    if (Y[ipt]<0) std::cout<<"!!!!! WEIGHT < 0 !!!!!"<<std::endl;
    xvals[index]=X[ipt];
    yvals[index]=Y[ipt];
    index++;
      //    } 
  }
  TString splineName = nameTag.Data();
  splineName.Append(Form("_spline_bin%d_par%d",ibin,isyst));
  checkSum--;
  if (!theSpline[ibin][isyst]) theSpline[ibin][isyst] = new TSpline3(splineName.Data(),xvals,yvals,index); 
  return;
}

hSplines::~hSplines()
{
}

hSplines::hSplines(TH1D* h, int nsyst, const char* name){
  if (name!=NULL) nameTag = name;
   nameTag = h->GetName();
   nHistoBins=h->GetNbinsX();
   nSyst=nsyst;
   //initialize base histogram
   baseHisto=h;
   cout<<"hSpline: Initialized splines for base histogram: "<<baseHisto->GetName()<<endl;
   //initialize modified histogram
   TString modhname = baseHisto->GetName();
   modhname.Append("_modified");
   int hbins = baseHisto->GetNbinsX();
   double hmin = baseHisto->GetBinLowEdge(1);
   double hmax = baseHisto->GetBinLowEdge(baseHisto->GetNbinsX())+baseHisto->GetBinWidth(1);
   modHisto = (TH1D*)baseHisto->Clone(modhname.Data());

   int basenbins = baseHisto->GetNbinsX();
   double basexmin  = baseHisto->GetBinLowEdge(1);
   double  basexmax = baseHisto->GetBinLowEdge(basenbins);
   checkSum = nHistoBins*nSyst;
}
