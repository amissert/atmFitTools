
#ifndef CALCENU_H
#define CALCENU_H

#include "TVector3.h"
#include "TMath.h"

using namespace std;

float calcEnu(float lpmom, TVector3 ldir, int muflg);


#ifdef CINTMODE
#include "calcEnu.cxx"
#endif


#endif
