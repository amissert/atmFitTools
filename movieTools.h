#ifndef  movieTools_H
#define  movieTools_H

#include <vector>
#include "TMath.h"
#include "TGraph.h"

using namespace std;

vector<double>* getParList(int nptstot, int nperiod, double xi, double delta, bool positive);
TGraph* plotV(const vector<double>* v);

#ifdef CINTMODE
#include "movieTools.cxx"
#endif

#endif  

