#ifndef MARKOVDIFF_H
#define MARKOVDIFF_H

#include "TRandom2.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include <math.h>
#include "TTree.h"
#include "TChain.h"
#include "shared.h"
#include <iostream>
#include "TFile.h"

//#define NMCMCPARS 500


/////////////////////////////////////////////////////
// Class to fill a TTree of MCMC step differences
class markovDiff{
  public:

  // constructor
  markovDiff(const char* fname, int burn=0,int chainflg=0);

  // internal vars
  TFile* mcmcfile;
  TFile* outfile;
  TTree* mcmcpars;
  TTree* diffpars;
  int nburn;
  float par[NMCMCPARS]; 
  float par1[NMCMCPARS]; //< vector of pars at point 1
  float par2[NMCMCPARS]; //< vector of pars at point 2
  float pardiff[NMCMCPARS]; //< difference vector between points
  int   parindex[NMCMCPARS]; //< atmFitPars index of each parameter
  int npars;
  TString outFileName; //< name of the output file
 
  // methods
  void fillDiffPars(int npairs); // fill the TTree of parameter difference vectors
  void fillDiffPars2(int npairs); // fill the TTree of parameter difference vectors
  void setaddresses();
  void setouttree();

};

#ifdef CINTMODE
#include "markovDiff.cxx"
#include "globalRandom.cxx"
#endif

#endif
