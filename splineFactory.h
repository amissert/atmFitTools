<<<<<<< HEAD
#ifndef SPLINEFACTORY_H
#define SPLINEFACTORY_H

#include "hSplines.h"

#include "shared.h"
#include "histoManager.h"
#include "sharedPars.h"
#include "atmFitPars.h"
#include "getSystWeight.cxx"


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//class for creating splines
class splineFactory{
  public:

  //constructors
  splineFactory(int nsamp, int nbin, int ncomp, int natt, int nsyst, const char* name, int nmode = 0, bool separateneutmode = false);  
  splineFactory(const char* parfile, bool separateneutmode = false);//< initialize using parameters in par file  
  splineFactory(){;}
  
  //internal variables
  TString parFileName; //< name of parameter file
  TString nameTag; //< set in constructor. this is the prefix for the output file
  TTree* mcTree; 
  fqProcessedEvent* mcEvt;
  histoManager* hManager; //manages all default histograms
  TH1D* hMC[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX][NPTSMAX]; //array for modified histograms for spline creation
  atmFitPars *fitPars;
  TH1D *hMCMode0[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX][NPTSMAX];
  TH1D *hMCMode1[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX][NPTSMAX]; 
  TH1D *hMCMode2[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX][NPTSMAX]; 
  TH1D *hMCMode3[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX][NPTSMAX]; 
  TH1D *hMCMode4[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX][NPTSMAX]; 
  TH1D *hMCMode5[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX][NPTSMAX]; 
  TH1D *hMCMode6[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX][NPTSMAX]; 
  TH1D *hMCMode7[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX][NPTSMAX]; 
  TH1D *hMCMode8[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX][NPTSMAX]; 
  TH1D *hMCMode9[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX][NPTSMAX]; 
  void resetModHistos();
  TFile *fout; //output file
  TString foutName;
  TH1D* htmp; //temporary histogram pointer
  int nSamp;
  int nBin;
  int nComp;
  int nAtt;
  int nSyst;
  TString sysParType; //< code denoting the type of parameterization used, see setupSysPar
  atmFitPars* fitPars;
  int nMode;
  bool separateNeutMode;
  double sysPar[NSYSTMAX]; //systematic parameter values
  double sysUnc[NSYSTMAX];  //systematic parameter uncertainties
  std::string sysName[NSYSTMAX];
  int nP[NSYSTMAX];
  double attribute[NATTMAX];
  double eventWeight;
  sharedPars* runpars; //< runtime parameters

  //for output tree
  TTree* splineTree;
  int nbin;
  int ncomponent;
  int nattribute;
  int nsample;
  int nmode;
  int nsystpar;
  int npoints;
  int nhistobins;
  double systParValues[NPTSMAX];
  double binWeight[NPTSMAX][NHBINSMAX];
  //methods
  double getEvtWeight(fqProcessedEvent* mcevent,int ipar,double value); //
  const static double sigvals[13];
  const static double maqevals[21];
  const static double binarys[2];
  double getEvtWeight(int ipar); //returns event weight after applying syst. par. 
  void setOutputFileName(const char* name){foutName=name;}
  TString getOutputFileName(){return foutName;}
  //
  void fillHistograms(int ipt,int isyst); //fills all histograms given weight
  void  makeManagerFromFile(const char* fname); //reads in histograms from histoFactory
  void fillBranches(int nsamp,int nbin,int ncomp,int natt,int isyst); //fills leaves of output tree
  void setMCTree(TTree* tr);
  void fillLeaves(int nsamp,int nbin,int ncomp,int natt,int isyst, int imode = -1); //fills leaves of output tree
  void setMCTree(TChain* tr);
  //build the splines
  void buildTheSplines();

  //debugging
  void debugtest();

  //do everything
  void runSplineFactory();

//  private:
  void fillAttributes();
  void incrementSystPars(int isyspar, double nsig);
  int getBest2RFitID();
  void setupHistos();
  void setupSystPars(); //sets up systematic parameters
  void incrementSystPars(double nsig);
  void incrementSystPars(double nsig, int i);
  void setAtmFitPars(atmFitPars *a);

};



#endif


#ifndef SPLINEFACTORY_C
#include splineFactory.cxx
#endif


