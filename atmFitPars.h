#ifndef __ATMFITPARS_H__
#define __ATMFITPARS_H__

#include "sharedPars.cxx"
#include "shared.h"
#include "TRandom2.h"
#include "sharedPars.h"
#include "covXsec.h"
#include "covBANFF.h"
#include "covBase.h"

using namespace std;

//class containing all atmospheric fit pars
class atmFitPars{
  public:
  
  ///////////////////////////////////////////////////////////////////
  //constructors
  atmFitPars(int isamp, int ibin, int icomp, int iatt, int nsyst=0);
  atmFitPars(int isamp, int ibin, int icomp, int iatt, const char* systype); 
  // use this constructor
 // atmFitPars(const char* parfile); //constructs from parameter file
  atmFitPars(const char* parfile, covBase *covm = 0); //constructs from parameter file

  ///////////////////////////////////////////////////////////////
  //numbers of various parametrs
  int nSamples;
  int nBins;
  int nComponents;
  int nAttributes;
  int nSysPars;
  int nTotPars;
  int flgUseNormPars;
  sharedPars* runpars;
  int nModes;
  sharedPars* runpars;
  covBase *cov;
  TRandom3 *rnd;

  ///////////////////////////////////////////////////////////////////
  //parameter values
  double histoNorm[NSAMPMAX][NBINMAX];
  double histoPar[NBINMAX][NCOMPMAX][NATTMAX][2];
  double histoParUncLo[NBINMAX][NCOMPMAX][NATTMAX][2];
  double histoParUncHi[NBINMAX][NCOMPMAX][NATTMAX][2];
  double sysPar[NSYSPARMAX];
  double sysParDefault[NSYSPARMAX];
  double sysParUnc[NSYSPARMAX];
  double pars[4000]; //< current values
  double parUnc[4000];
  int   fixPar[4000]; //< array of fix flags for parameters
  double bestpars[4000];
  int   parIndex[NBINMAX][NCOMPMAX][NATTMAX][2]; //< stores 1D array position for bias/smear pars
  int   sysParIndex[NSYSPARMAX]; //< stores 1D array position for systematic pars
  int   normParIndex[NSAMPMAX][NBINMAX]; //< stores 1D array position for normalization pars
  double sysParNom[NSYSPARMAX];
  double sysParUnc[NSYSPARMAX];
  double sysParUp[NSYSPARMAX];
  double sysParLow[NSYSPARMAX];
  string sysParName[NSYSPARMAX];
  double norm;  
  double parsProp[4000];
  double bestpars[4000];
  int   parIndex[NBINMAX][NCOMPMAX][NATTMAX][2];
  float fScale;

  //////////////////////////////////////////////////////////////
  //methods
  void setNorm(double x){norm=x;}
  void initPars(const char* systype=""); //< sets parameters to initial values
  int getParIndex(int ibin, int icomp, int iatt, int imod){return parIndex[ibin][icomp][iatt][imod];}
  double getParameter(int ipar){return pars[ipar];}
  double getHistoParameter(int ibin, int icomp, int iatt, int imod);
  double getSysParameter(int isys);
  double getNormParameter(int isamp, int ibin);
  void setParameter(int ipar, double value);
  void setSysParameter(int ipar, double value);
  void setParameter(int ibin, int icomp, int iatt, int imod, double value); 
  void setSysParUnc(int isys,double value){sysParUnc[isys]=value;}

  //////////////////////////////////////////////////////////////
  //methods
  void proposeStep();
  void acceptStep();
  void setStepSize(float f) {fScale = f; cov->setStepScales(f);}
  void setSeed(int i) {rnd->SetSeed(i);}
  double getParameter(int ipar){return pars[ipar];}
  double getPropParameter(int ipar) { return parsProp[ipar]; }
  double* getParameters() {return pars;}
  void fixParameter(int ipar);
  void fixParameter(int ibin,int icomp,int iatt, int imod);
  void fixAllSmearPars(int isfixed=1);
  void setRandSysPar(); //sets systematic parameters to random values
  int  checkFixFlg(int ibin,int icomp,int iatt, int imod);
  void resetDefaults();
  void printParValues();
  void setCov(covBase *covariance);

  int binOfPar[4000];
  int compOfPar[4000];
  int attOfPar[4000];
  int typeOfPar[4000];
  std::string sysType;

  //saving and reading pars
  void savePars(const char* filename);
  void readPars(const char* filename);
  void printPars();
};

#endif

#ifndef ATMFITPARS_C
#include "atmFitPars.cxx"
#endif


