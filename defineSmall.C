#ifndef DEFINESMALL_CXX
#define DEFINESMALL_CXX


#include "TMath.h"
#include "TH1D.h"
#include "TRandom2.h"
#include "stats.C"



//////////////////////////////////////////////////////////////
// Define a "small" variation for standard FOM using MC
double defineSmallSignificance(double Ns, double Nb, int npts){

  // make random generator
  TRandom2* randgen = new TRandom2(npts);

  // array fo values
  const int N = npts;
  double signif[N];

  // fill array
  for (int i=0; i<npts; i++){

    // get random signal
    double ss = randgen->PoissonD(Ns);
//    cout<<"ss: "<<ss<<endl; 

    // get random bg
    double bg = randgen->PoissonD(Nb);
//    cout<<"bg: "<<bg<<endl; 
    
    double sig = 0.;
    if (ss+bg > 0.) sig = (ss*ss)/(ss+bg);
//    cout<<"sig: "<<sig<<endl; 

    // fill array
    signif[i]  = sig;
    
  }

  // get variance 
  double var = arrayvarD(signif,npts);
  return TMath::Sqrt(var);

}

 
/*
//////////////////////////////////////////////////////////////
// Define a "small" variation for FOM using MC in a single bin
double defineSmall(double Ns, double Nb, double avgpow, double avgsigsys, double avgbgsys, int npts=10000){

  // make random generator
  TRandom2* randgen = new TRandom2(npts);

  // array fo values
  const int N = npts;
  double signif[N];

  // fill array
  for (int i=0; i<npts; i++){

    // get random signal
    double ss = randgen->PoissonD(Ns);
    double pow = avgpow*ss;

    // get random bg
    double bg = randgen->PoissonD(Nb);
   
    // get systematic component
    double sys_ss = avgsigsys*ss;
    double sys_bg = avgbgsys*bg;

    double numer = (pow*pow);
    double denom = (ss + bg + ((sys_ss + sys_bg)*(sys_ss*sys_bg)));    

    double fom = 0.;
    if (denom > 0.) fom = numer/denom;

    // fill array
    signif[i]  = fom;
    
  }

  // get variance 
  float var = arrayvarD(signif,npts);
  return TMath::Sqrt(var);

}
*/


float getVariation(float pow, float nev, float sys, float powavg, float sysavg){

  float fom0 = (pow*pow)/(nev + (sys*sys));
//  float fom1 = ((pow-1.)*(pow-1.))/((nev-1.) + ((sys-0.5)*(sys-0.5)));
//  float fom2 = ((pow+1.)*(pow+1.))/((nev+1.) + ((sys+0.5)*(sys+0.5)));
//  float powavg = 0.5;
//  float fom1 = ((pow-powavg)*(pow-powavg))/((nev-1) + ((sys-sysavg)*(sys-sysavg)));
//  float fom2 = ((pow+powavg)*(pow+powavg))/((nev+1) + ((sys+sysavg)*(sys+sysavg)));

  float fom1 = ((pow-powavg)*(pow-powavg))/((nev) + ((sys)*(sys)));
  float fom2 = ((pow+powavg)*(pow+powavg))/((nev) + ((sys)*(sys)));

  return TMath::Abs((fom2-fom1))/2.;
}


float getSysVariation(float pow, float nev, float sys){

  float scaleup = 1.025;
  float scaledown = 0.975;
  float fom1 = ((pow)*(pow))/((nev) + ((sys*scaleup)*(sys*scaleup)));
  float fom2 = ((pow)*(pow))/((nev) + ((sys*scaledown)*(sys*scaledown)));

  return TMath::Abs((fom2-fom1))/1.;
}


///////////////////////////////////////////////////////////////////////
// one way to define small
float defineSmall(int nbins, float* pow, float* nev, float* sys){

  // first find maximum of figure of merit
  float fomax = 0.;
  int maxbin = -1;
  for (int i=0; i<nbins; i++){
    float fom = 0.;
    if ((nev[i]+sys[i])>0.){
      fom = (pow[i]*pow[i])/(nev[i] + (sys[i]*sys[i]) );
    }
    if (fom>fomax){
      fomax = fom;
      maxbin = i;
    }
  }
  cout<<"maxbin: "<<maxbin<<endl;
  cout<<"fomax: "<<fomax<<endl;

  // now, get average power estimate in that bin
  float powavg = TMath::Abs(pow[maxbin])/nev[maxbin];
  float sysavg = TMath::Abs(sys[maxbin])/nev[maxbin];

  cout<<"powavg: "<<powavg<<endl;

  // see how much a variation changes
  float var = getVariation(pow[maxbin],nev[maxbin],sys[maxbin],powavg,sysavg);

  return var;
  
}

///////////////////////////////////////////////////////////////////////
// anothe rway to define small
float defineSmall2(int nbins, float* pow, float* nev, float* sys){

  // first find maximum of figure of merit
  float fomax = 0.;
  int maxbin = -1;
  for (int i=0; i<nbins; i++){
    float fom = 0.;
    if ((nev[i]+sys[i])>0.){
      fom = (pow[i]*pow[i])/(nev[i] + (sys[i]*sys[i]) );
    }
    if (fom>fomax){
      fomax = fom;
      maxbin = i;
    }
  }

  // see how much a variation changes
  float var = getSysVariation(pow[maxbin],nev[maxbin],sys[maxbin]);

  return var;
  
}


double getAbsDifference(TH1D* h1, TH1D* h2){

  int nbins = h1->GetNbinsX();
  double diffsum = 0.;

  for (int i=0; i<=nbins; i++){

//    double absdiff = TMath::Abs(h1->GetBinContent(i) - h2->GetBinContent(i) );
    double absdiff = (h2->GetBinContent(i) - h1->GetBinContent(i) );
    diffsum += absdiff;

  }

  return diffsum;

}










#endif
