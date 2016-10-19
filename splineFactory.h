#ifndef SPLINEFACTORY_H
#define SPLINEFACTORY_H

#include "hSplines.h"

#include "shared.h"

#include "histoManager.h"
#include "sharedPars.h"
#include "atmFitPars.h"
#include "getSystWeight.h"

#define NHBINSMAX 300
#define NPTSMAX   5



/////////////////////////////////////////////////////////////////////////////////////////////////////////
//class for creating splines
class splineFactory{
  public:

  //constructors
  splineFactory(int nsamp, int nbin, int ncomp, int natt, int nsyst, const char* name);
  splineFactory(const char* parfile);//< initialize using parameters in par file
  splineFactory(){;}
  
  //internal variables
  TString parFileName; //< name of parameter file
  TString nameTag; //< set in constructor. this is the prefix for the output file
  TTree* mcTree; 
  fqProcessedEvent* mcEvt;
  histoManager* hManager; //manages all default histograms
  TH1D* hMC[NSAMPMAX][NBINMAX][NCOMPMAX][NATTMAX][NPTSMAX]; //array for modified histograms for spline creation
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
//  double sysPar[NSYSTMAX]; //systematic parameter values
//  double sysUnc[NSYSTMAX];  //systematic parameter uncertainties
  atmFitPars* fitPars;
  double attribute[NATTMAX];
  double eventWeight;
  sharedPars* runpars; //< runtime parameters
  double sigmaValues[NPTSMAX];

  //for output tree
  TTree* splineTree;
  int nbin;
  int ncomponent;
  int nattribute;
  int nsample;
  int nsystpar;
  int npoints;
  int nhistobins;
  double systParValues[NPTSMAX];
  double binWeight[NPTSMAX][NHBINSMAX];
  TH2D* h2DWeights;

  //methods
  //this needs to be modified for each systematic paramater to add
 // double getEvtWeight(int ipar); //returns event weight after applying syst. par. 
  double getEvtWeight(fqProcessedEvent* mcevent,int ipar,double value); //
  void setOutputFileName(const char* name){foutName=name;}
  TString getOutputFileName(){return foutName;}
  //
  void  makeManagerFromFile(const char* fname); //reads in histograms from histoFactory
  void fillBranches(int nsamp,int nbin,int ncomp,int natt,int isyst); //fills leaves of output tree
  void setMCTree(TTree* tr);
  void setupOutputFile();

  //build the splines
  void buildTheSplines();
  void buildSplineForPar(int isyspar);
  //debugging
  void debugtest();

  //do everything
  void runSplineFactory();

  void setupHistos();
//  void setupSystPars(); //sets up systematic parameters
  void incrementSystPars(int isyspar, double nsig);
};


#ifdef CINTMODE
#include "splineFactory.cxx"
#endif


#endif



