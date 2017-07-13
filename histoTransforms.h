#ifndef HISTOTRANSFORMS_H
#define HISTOTRANSFORMS_H

#include "TH1D.h"
#include "TRandom2.h"
#include "TMath.h"
#include "time.h"
#include <iostream>
#include "TCanvas.h"
#include "TGraph.h"
#include "shared.h"

//#ifndef GLOBAL_RANDOM
//#define GLOBAL_RANDOM
//TRandom2* randy = new TRandom2();
//#endif

void smearThisGraph(TGraph* gr, double smear, double bias);
void smoothGraph(TGraph* gr);
double gIntegral(TGraph* gr, double xmin, double xmax, int sampling=25);
void shiftGraph(TGraph* gr, double smear, double bias);
TH1D* testBumpD(int nev,double sig=1.0,double mean=0.0,const char* name="testbumb");
double graph2histo(TGraph* gr, TH1D* h);
void applyLoBound(TGraph* gr, TH1D* h, double lobound);
TGraph* histo2graph(TH1D*h);
double B(double x,double a, double b);
double getNoiseFactor(TH1D* hh);
void mySmooth2(TH1D* hh,double factor=3.0);
void mySmooth(TH1D* hh,double factor=3.0);
void rebinHisto(TH1D* holdbin, TH1D* hnewbin);
TH1D* testBump(int nev,double sig=1.0,double mean=0.0,const char* name="testbumb");
TH1D* testTable(int nev,double width=4.0,double mean=0.0);
void smearThisHisto(TH1D &hh, double spread, double bias=0.);
void compareEvtByEvt(int nevts, double mean0, double sig0,double scale, double bias);
void smearThisHistoFast(TH1D &hh, double* hcontent, double spread,  double bias, double normscale=1.);
void smearThisHistoFastBias(TH1D &hh, double* hcontent, double bias, double normscale=1.);
void smearThisHistoFastMean(TH1D &hh, double* hcontent, double spread, double mean, double bias, double normscale=1.);
double testtime();

#ifdef CINTMODE
#include "histoTransforms.cxx"
#endif

#endif
