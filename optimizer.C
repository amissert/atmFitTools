#include "TH2D.h"
#include "TH1D.h"
#include "TMath.h"
#include "fQreader.C"
#include "TTree.h"
#include "TFile.h"
#include "classify.C"
#include "TH2FV.C"
#include <iostream>
#include "TGraph2D.h"
#include "TGraph.h"


using namespace std;

class optimizer{
  public:

  // constructor
  optimizer();

  // internal vars
  TH2FV* hsys;
  TTree* tree;
  TH2D* hfom;
  TH2FV* hsyscopy;
  int nevents;
  classify* classifier;
  fQreader* fqevent;
  TGraph2D* gmap;
  TGraph*   gbest;
  TH2FV*    hmap;
  double normfactor;
  void setsyshisto(TH2D* hh);
  void setsyshisto(const char* fname, const char* hname);
  double calcfom(double towall, double wall);
  double getsyserror(double towall, double wall);
  double getsyserror(fQreader* fqevent);
  void init();
  void mapfom(int npts, double xmin, double xmax, double ymin, double ymax);
  // methods
  
  
};

/////////////////////////////////
// map out the fom in X (towall) and Y (wall)
void optimizer::mapfom(int npts, double xmin, double xmax, double ymin, double ymax){

  // increments towall and wall
  double dx = (xmax-xmin)/(double)npts;
  double dy = (ymax-ymin)/(double)npts;

  // set initial values
  double xx;
  double yy;
  double zz;

  // arrays for TGraph2D
  const int N = npts*npts;
  double X[N];
  double Y[N];
  double Z[N];
  int nn = 0; //< counts total number of points;
  double bestvalue = -1.; //< best (max) value of f.o.m.
  double bestwall = 0.;
  double besttowall = 0.;

  // set up hmap
  hmap = new TH2FV("hmap",-1,npts,xmin,xmax,npts,ymin,ymax);

  // loop over points and calculate f.o.m.
  xx = xmin+0.1;
  for (int ibin=1; ibin<hmap->GetNumberOfBins(); ibin++){

     // get (towall, wall) cut values
     xx = hmap->GetBinCenterX(ibin);
     yy = hmap->GetBinCenterY(ibin);
 
     // skip nonsense
     if (yy>xx) continue;
     
     // calculate figure of merit for test cut
     zz = calcfom(xx,yy);

     if (zz>bestvalue){
      bestvalue = zz;
      bestwall = yy;
      besttowall = xx;
    }

    hmap->SetBinContent(ibin,zz);
    
  }
  double xbest[2];
  double ybest[2];
  xbest[0] = besttowall;
  ybest[0] = bestwall;
  gbest = new TGraph(1,xbest,ybest);
  hmap->SetContour(100);
  hmap->Draw("colz");
  gbest->Draw("same*");

  return;
}


/*
// mmap out the fom
void optimizer::mapfom(int npts, double xmin, double xmax, double ymin, double ymax){
  // increments towall and wall
  double dx = (xmax-xmin)/(double)npts;
  double dy = (ymax-ymin)/(double)npts;
  // set initial values
  double xx;
  double yy;
  // arrays for TGraph2D
  const int N = npts*npts;
  double X[N];
  double Y[N];
  double Z[N];
  int nn = 0; //< counts total number of points;
  double bestvalue = -1.; //< best (max) value of f.o.m.
  double bestwall = 0.;
  double besttowall = 0.;
  // set up hmap
  hmap = new TH2FV("hmap",-1,npts,xmin,xmax,npts,ymin,ymax);
  // loop over points and calculate f.o.m.
  xx = xmin+0.1;
  for (int ipt=0; ipt<npts; ipt++){
    yy = ymin+0.01;
    for (int jpt=0; jpt<npts; jpt++){
      if (yy>xx) continue;
      X[nn] = xx;
      Y[nn] = yy;
      Z[nn] = calcfom(xx,yy);
      if (Z[nn]>bestvalue){
        bestvalue = Z[nn];
        //bestbin = nn;
        bestwall = Y[nn];
        besttowall = X[nn];
      }
      // hmape
      int ibin = hmap->Fill(xx,yy);
      hmap->SetBinContent(ibin,Z[nn]);
      yy += dy;
      nn++;
    }
    xx += dx;
  }
  // make TGraph
  gmap = new TGraph2D(nn,X,Y,Z);
//  gmap->SetContour(100);
 // gmap->Draw("colz");
 // best points
  double xbest[2];
  double ybest[2];
  xbest[0] = besttowall;
  ybest[0] = bestwall;
  gbest = new TGraph(1,xbest,ybest);
  hmap->SetContour(100);
  hmap->Draw("colz");
  gbest->Draw("same*");
}
*/


// set a wall v towall histogram to calculate the systematic errors
void optimizer::setsyshisto(const char* fname, const char* hname){
  TFile* file = new TFile(fname);
  hsys = (TH2FV*)file->Get(hname);
  cout<<"optimizer: Setting systematics histogram to "<<hname<<" from "<<hname<<endl;
  hsyscopy = (TH2FV*)hsys->Clone("hsyscopy");
  hsyscopy->Reset("ICE");
  return;
}

// get the systematic error for this event
// "slow" method since we are using hsys to get bin
// should just do this once in preprocessing
double optimizer::getsyserror(double towall, double wall){
  //int fvbin = fqevent->nbin+1;
  int fvbin = hsyscopy->Fill(towall,wall);
  double syserr = hsys->GetBinContent(fvbin);
  return syserr;
}

double optimizer::getsyserror(fQreader* fqevent){
  return getsyserror(fqevent->fq1rtowall[0][2],fqevent->fq1rwall[0][2]);
}

void optimizer::init(){
  classifier = new classify(0);
  fqevent = new fQreader(tree);
  hsyscopy = (TH2FV*)hsys->Clone("hsyscopy");
  classifier = new classify(); 
  hsyscopy->Reset("ICES");
  nevents = tree->GetEntries();
  return;
}
  

double optimizer::calcfom(double towall, double wall){
  
  // loop over events and calc fom
  double Ns = 0;
  double Nb = 0;
  double Sys = 0;
  for (int ievt=0; ievt<nevents; ievt++){
   
    if ((ievt%5000)==0){
      int pct = (ievt*100)/nevents;
      cout<<" ..."<<pct<<"%"<<" done"<<endl;
    }
    // read in event data
    tree->GetEntry(ievt);

    // see if event passes FV cuts for this point
    if (fqevent->fq1rwall[0][1]<wall) continue;
    if (fqevent->fq1rtowall[0][2]<towall) continue;
    if (classifier->ismuonclass(fqevent)<0) continue;
    // see if event passes muon cuts 
    if (TMath::Abs(fqevent->mode)==1){
      Ns+=normfactor;
    }
    else{
      Nb+=normfactor;
    }

    // add up the systematic error contribution
    double eventsys = getsyserror(fqevent);
    Sys+=(eventsys*normfactor);
  }

  // get the total number of events in your sample
  double Ntot = Ns+Nb;

  // debugging
  cout<<"N: "<<Ntot<<endl;
  cout<<"Ns: "<<Ns<<endl;
  cout<<"Nb: "<<Nb<<endl;
  cout<<"Sys: "<<Sys<<endl;

  // return f.o.m.
  double fom = Ns/TMath::Sqrt(Ntot+(Sys*Sys));

  return fom;
}

optimizer::optimizer(){
  normfactor = 1.0;
  classifier = new classify(); 
}








