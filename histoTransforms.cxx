#ifndef HISTOTRANSFORMS_C
#define HISTOTRANSFORMS_C

#include "histoTransforms.h"

using namespace std;

////////////////////////////////////////////////////
// Apply smearing/bias to graph
void smearThisGraph(TGraph* gr, double smear, double bias){
 

  double* Y = gr->GetY();
  double* X = gr->GetX();
  int npts = gr->GetN();

  for (int ipt=0; ipt<npts; ipt++){
    gr->SetPoint(ipt, (smear*X[ipt] + bias), Y[ipt]); 
  }
  
  
  return;
}

////////////////////////////////////////////////////
// Smooth over a graph
void smoothGraph(TGraph* gr){

  // smooth weights
  double ww0 = 1.0;
  double ww1 = 6.06530659712633424e-01;
  double ww2 = 1.35335283236612702e-01;
  double wsum = ww0 + (2*ww1) + (2.*ww2);
  double* X = gr->GetX();
  double* Y = gr->GetY();
  int N = gr->GetN();

  // first point remians same
  // second point is avg
  double yy = (ww0*Y[1] + ww1*Y[0] + ww1*Y[2])/(ww0+ww1+ww1);
  gr->SetPoint(1,X[1],yy);
  // begin loop on third point
  for (int ipoint=2; ipoint<=(N-2); ipoint++){
    yy = ((ww0*Y[ipoint]) + (ww1*Y[ipoint-1]) + (ww1*Y[ipoint+1]) + (ww2*Y[ipoint-2]) + (ww2*Y[ipoint+2]) )
         / wsum; 
    gr->SetPoint(ipoint,X[ipoint],yy);
  }
  // avg. 2nd last point
  yy = (ww0*Y[N-1] + ww1*Y[N-2] + ww1*Y[N])/(ww0+ww1+ww1);
  gr->SetPoint(1,X[1],yy);
  
  return;
}

////////////////////////////////////////////////////
// Integrate a graph
double gIntegral(TGraph* gr, double xmin, double xmax, int sampling){

  // step size
  double dx = (xmax-xmin)/((double)sampling);

  // add areas
  double area = 0.;
  double xx   = xmin;
 
  // get X boundaries
  double* X = gr->GetX();
  double loBound = X[0];
  double hiBound = X[gr->GetN()-1];

  for (int i=0; i<sampling; i++){
    // rectangles
    double xvalue = xx + (dx/2.);
    // make sure in bounds
    if ((xvalue<loBound)||(xvalue>hiBound)) continue;
    //rectangles
    //area += dx*gr->Eval( xvalue,0,"s" ); 
    //area += dx*gr->Eval( xvalue ); 
    area += dx*gr->Eval( xvalue ); 
    xx+=dx;
  }

  return area;

}

void shiftGraph(TGraph* gr, double smear, double bias){
  int N = gr->GetN();
  double* X = gr->GetX();
  double* Y = gr->GetY();
  for (int i=0; i<N; i++){
     gr->SetPoint(i,smear*X[i] + bias,Y[i]/smear); 
  }

  return;
}

////////////////////////////////////////////////////
// converts graph to histogram
// histogram bin contents will be re-written
double graph2histo(TGraph* gr, TH1D* h){
  
  // clear bin contents
  h->Reset();

  // calc loss
  double *X = gr->GetX();
  double binw = h->GetBinWidth(1);
  int    nbins = h->GetNbinsX();
  double hxmin = h->GetBinLowEdge(0);
  double hxmax = h->GetBinLowEdge(nbins) + binw;

  // fill bin contents from graph
  for (int ibin=1; ibin<=nbins; ibin++){
    double xmin = h->GetBinLowEdge(ibin);
    double xmax = xmin + binw;
    double area = gIntegral(gr,xmin,xmax);
    h->SetBinContent(ibin,area);
    h->SetBinError(ibin,TMath::Sqrt(area));
  }
 
  return 1.;
}


/////////////////////////////////////////////////////////////
// Modify histogram filled by graph with physical lower bound
// Can be used as a template for enforcing an upper bound as well
void applyLoBound(TGraph* gr, TH1D* h, double lobound){
  
  // get bin width (assume constant);
  double binw = h->GetBinWidth(1);
  // identify the critical bin where the lower bound lies
  int critbin = -1;
  for (int ibin=1; ibin<=h->GetNbinsX(); ibin++){

     double binloval = h->GetBinLowEdge(ibin);
     double binhival  = binloval + h->GetBinWidth(ibin); 
     if ((binloval <= lobound) && (binhival> lobound)){
       critbin = ibin;
       break;
     }

  }

  // integrate un-physical bins
  double binsum=0;
  for (int ibin=(critbin-1); ibin>=0; ibin--){
    binsum+=h->GetBinContent(ibin);
  }
 
  // add to critical bin
  h->SetBinContent(critbin,h->GetBinContent(critbin)+binsum);

  // clear un-physical bins
  for (int ibin=(critbin-1); ibin>=0; ibin--){
    h->SetBinContent(ibin,0.);
  }
  
  //
  return;
}




/////////////////////////////////////////////////////
// Similar to TGraph constructor with with more care
// at histogram endpoints
TGraph* histo2graph(TH1D*h){
  
  // graph paramters
  int nbins =  h->GetNbinsX();
  double binw = h->GetBinWidth(1);
  const int N = nbins + 2;
  double X[N];
  double Y[N];

  // set first point
  X[0] = h->GetBinLowEdge(1);
  Y[0] = 0.;

  // set last point
  X[N-1] = h->GetBinLowEdge(nbins)+binw;
  Y[N-1] = 0.;
  
  // set other points
  for (int i=1; i<=nbins; i++){
    X[i] = h->GetBinCenter(i);
    Y[i] = h->GetBinContent(i);
  }

  TGraph* g = new TGraph(N,X,Y);

  return g;
}

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




///////////////////////////////////////////////////
// Rebin using intermediate TGraph
void rebinHisto(TH1D* holdbin, TH1D* hnewbin){
  
  // make tgraph using old bins
  TGraph* gr = histo2graph(holdbin);

  // integrate into new bins
  graph2histo(gr,hnewbin);

  gr->Delete();

  // done
  return;
}




///////////////////////////////////////////////////////////////////////////////////
//Custom smoothing method
void mySmooth2(TH1D* hh,double factor){

  //////////////////////////////
  //get adjecent bin weights
  
  // use Gaussian weights
  double sigma = factor*getNoiseFactor(hh); //< set sigma using noise
  if (sigma==0) return; //no events in histogram
  sigma = 1.0;
//  double w2 = TMath::Gaus(2,0,sigma,1);
//  double w1 = TMath::Gaus(1,0,sigma,1);
//  double w0 = TMath::Gaus(0,0,sigma,1);

  double w2 = TMath::Gaus(4,0,sigma,1);
  double w1 = TMath::Gaus(2,0,sigma,1);
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



///////////////////////////////////////////////////////////////////////////////////
//Custom smoothing method
void mySmooth(TH1D* hh,double factor){

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
TH1D* testBumpD(int nev,double sig,double mean,const char* name){
  TH1D* hset;
  hset->SetDefaultSumw2(kTRUE);
  TH1D* h = new TH1D(name,name,40,-10,10);
  for (int  i=0;i<nev;i++){
    h->Fill(randy->Gaus(mean,sig));
  }
  return h;
}

TH1D* testBump(int nev,double sig,double mean,const char* name){
  TH1D* hset;
  hset->SetDefaultSumw2(kTRUE);
  TH1D* h = new TH1D(name,name,40,-10,10);
  for (int  i=0;i<nev;i++){
    h->Fill(randy->Gaus(mean,sig));
  }
  return h;
}

TH1D* testTable(int nev,double width,double mean){
  TH1D* h = new TH1D("testbump","testbump",50,-10,10);
  for (int  i=0;i<nev;i++){
    h->Fill((randy->Rndm()*width) - (width/2.)+mean);
  }
  return h;
}



//////////////////////////////////////////////////////////////
//performs histogram modifications
//(old method, use TGraph methods now)
void smearThisHisto(TH1D &hh, double spread, double bias){

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
    hh.SetBinContent(newbin,sum);
    hh.SetBinError(newbin,TMath::Sqrt(sum));
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


//smear it faster without mean
void smearThisHistoFast(TH1D &hh, double* hcontent, double spread,  double bias, double normscale){

  int nbins=hh.GetNbinsX();
 // double oldintegral = hh.Integral();
  double binw = hh.GetBinWidth(1);

  //parameters for calculations
  double binedge;
  double sum;
  double smear = 1./spread;
  double binerr;
  double weight;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double sumw;

  //loop over bins and modify contents
  for (int newbin=1;newbin<=nbins;newbin++){
    sum = 0.;
    sumw = 0;
    binerr=0.;
    binedge = hh.GetBinLowEdge(newbin);
//    ymin = binedge - bias;
//    ymax = ymin + binw;
    ymin = ((binedge-bias)*smear);
    ymax = ymin + (binw*smear);
    for (int oldbin=1;oldbin<=nbins;oldbin++){
      xmin = hh.GetBinLowEdge(oldbin);
      xmax = (xmin+binw);
      weight =  B(xmax,ymin,ymax)-B(xmin,ymin,ymax);
      double binc = hcontent[oldbin];
      sum+=(weight*binc);
      binerr += weight*weight*binc;
      sumw += weight;
    }
    hh.SetBinContent(newbin,sum);
    hh.SetBinError(newbin,TMath::Sqrt(binerr));
  }
  return;
}



//smear it faster without mean
void smearThisHistoFastBias(TH1D &hh, double* hcontent, double bias, double normscale){

  int nbins=hh.GetNbinsX();
 // double oldintegral = hh.Integral();
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
  double sumw;

  //loop over bins and modify contents
  for (int newbin=1;newbin<=nbins;newbin++){
    sum = 0.;
    sumw = 0;
    binerr=0.;
    binedge = hh.GetBinLowEdge(newbin);
    ymin = binedge - bias;
    ymax = ymin + binw;
    for (int oldbin=1;oldbin<=nbins;oldbin++){
      xmin = hh.GetBinLowEdge(oldbin);
      xmax = (xmin+binw);
      weight =  B(xmax,ymin,ymax)-B(xmin,ymin,ymax);
      double binc = hcontent[oldbin];
      sum+=(weight*binc);
      binerr += weight*weight*binc;
      sumw += weight;
    }
    hh.SetBinContent(newbin,sum);
    hh.SetBinError(newbin,TMath::Sqrt(binerr));
  }
  return;
}



//smear it faster
void smearThisHistoFastMean(TH1D &hh, double* hcontent, double spread, double mean, double bias, double normscale){

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
  double binwidth;
  double sum;
  double binerr;
  double weight;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double smear = 1./spread;
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

