#ifndef HISTCOMPARE_H
#define HISTCOMPARE_H

#include "histoManager.C"
#include "TMath.h"
#include "TRandom2.h"
#include "histoTransforms.C"
#include "TFitter.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
using namespace std;


//class to compare histograms and evaluate likelihoods
class histoCompare{
  public:

  //constructors//
  histoCompare(const char* thename);  //standard constructor

  //internal variables
  TString nameTag;  //name associated with this instance
  int nSamp;  //number of samples
  int nBin;  //number of bins
  int nComp;  //number of  components
  int nAtt;  //nummboer of attributes
  float tunePar; //tuning parameter for MCMC
  //tools for histogram manager management
  //created histo manager from file
  void readFromFile(const char* rootname,int isamp,int ibin, int icomp, int natt);
  histoManager* hManager;

  //atmospheric pars
  atmFitPars* thePars;

  
  float Norm;
  float Par[NBINMAX][NCOMPMAX][NATTMAX][2];
  float sysPar[NSYSPARMAX];
  float sysParUnc[NSYSPARMAX];
  
  float fixPar[NBINMAX][NCOMPMAX][NATTMAX][2];
  float bestPar[NBINMAX][NCOMPMAX][NATTMAX][2];
 // float errParLo[NBINMAX][NCOMPMAX][NATTMAX][2];
 // float errParHi[NBINMAX][NCOMPMAX][NATTMAX][2];
  float errParLo[1000];
  float errParHi[1000];

  TString parName[NBINMAX][NCOMPMAX][NATTMAX][2];
  TString binName[NBINMAX];
  TString compName[NCOMPMAX];
  TString attName[NCOMPMAX];
  void setBinName(int ibin, const char* name){binName[ibin]=name;}
  void setCompName(int icomp, const char* name){compName[icomp]=name;}
  void setAttName(int iatt, const char* name){attName[iatt]=name;}
  void setupPars(int nsyspars=0); //sets up all parameters
  //post-fit toolts
  void profileL(int ibin, int icomp, int iatt, int imod, float range, int npts=100);
  void profileL(int ipar,float range, int npts=100);
  void showFitHisto(int isamp,int ibin,int icomp,int iatt);
  void showFitEffect(int isamp,int ibin,int icomp,int iatt);
  void showFitResult(int isamp,int ibin,int iatt);
  void showFitPars(int ibin,int iatt,int imod);
  void showModHiso(int isamp,int ibin, int icomp, int iatt, float smear, float bias);
  void runMCMC(int nsteps);
 // float getErrLo(int ibin,int icomp,int iatt,int imod);
  float getErrLo(int isyst);
  float getErrHi(int isyst);
  TH1F* getModifiedHisto(int ibin, int icomp, int iatt){return hManager->getSumHistogramMod(ibin,icomp,iatt);}
 // TH1F* hMod[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX];

  void setupSplines(const char* fname,int npar){hManager->readSplinesFromFile(fname,npar);}
  //tools for adding histogram directly..for debugging
  void addHistogram(TH1F* h,int dataflg);
  int  rebinFactor;
  void setRebinFactor(int ifact){rebinFactor=ifact;}
  TH1F* hData[10];
  TH1F* hMC[10];
  TH1F* hModDebug;
  TH1F* hMod;
  TH1F* hTmp; //temporary histogram container
  TH1F* hPar; //for fit parameters
  TGraphAsymmErrors* gPar;
  TH1F* hParErrLo;
  TH1F* hParErrHi;
  TH1F* hProf;
  TH1F* showSmear(TH1F* h, float smear, float bias);
  void showMod(int imchist);
  TH1F* hTot;
  TCanvas* cc;
  int nDataHist;
  int nMCHist;
  float parDebug[10][2];
//  float getLDebug(); 
  
  //likelihood evaluateions
  float getSumSq(TH1F* h1, TH1F* h2);
  float getLnL(TH1F* h1, TH1F* h2);
  float getNDiff();
  float getTotSumSq();
  float getTotLnL();
  static void sumSqWrapper(int& ndim, double* gout, double& result, double par[], int flg);
  static void lnLWrapper(int& ndim, double* gout, double& result, double par[], int flg);
  void  getTotLnL1D(float& result,int npar, float par[]);
  //for debuggint and play
  float getTotSumSqDebug();
  static void sumSqDebugWrapper(int& ndim, double* gout, double& result, double par[], int flg);
  static void nDiffDebugWrapper(int& ndim, double* gout, double& result, double par[], int flg);
  void minSumSqDebug();
  void minSumSq();
  void sumSqPrefit();
  void LnLFit();
  void LnLPreFit();
  void singleParFit(int ipar);
  void sysParFit();
  void drawResult(int ihist);
  void timetest(int ntry);
  void saveFitPars(const char* filename); //< write parameters and outputs to a file
  void readFitPars(const char* filename); //< read parameters from a file
  void tuneMCMC(int ncyles=1,int nsteps=150,float goal=0.25);
  //staticthis for fits
  static histoCompare* staticthis;
};

#endif
