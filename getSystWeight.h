#ifndef GETSYSTWEIGHT_H
#define GETSYSTWEIGHT_H

#include "atmFitPars.h"
#include "fqProcessedEvent.h"
#include "TString.h"

using namespace std;

double getSystWeight(const char* systype, fqProcessedEvent* mcevent, int ipar, double value);

#ifdef CINTMODE
#include "getSystWeight.cxx"
#endif
#endif
