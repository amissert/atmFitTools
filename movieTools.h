#ifndef  movieTools_INC
#define  movieTools_INC

#include <vector>
#include "TMath.h"
#include "TGraph.h"

using namespace std;


vector<double>* getParList(int nptstot, int nperiod, double xi, double delta, bool positive){

  double factor = 1.;
  if (!positive) factor = -1.;
  vector<double>* vptr = new vector<double>;
  for (int i=0; i<nptstot; i++){
    double phase = (double)i/(double)nperiod;
    double value = xi + factor*(delta*TMath::Cos(2.*TMath::Pi()*phase));
    vptr->push_back(value);     
  }

  return vptr;

}


TGraph* plotV(const vector<double>* v){

  const int N = v->size();
  double X[N]; 
  double Y[N]; 
  for (int i=0; i<N; i++){
    X[i] = (double)i;
    Y[i] = v->at(i);
  }
  TGraph* g = new TGraph(N,X,Y); 

  g->Draw("alp");

}

#endif  

