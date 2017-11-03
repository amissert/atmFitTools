#ifndef MARKOV_TOOLS_H
#define MARKOV_TOOLS_H

#include "TRandom2.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include <math.h>
#include "TTree.h"
#include <iostream>
#include "TFile.h"
#include "atmFitPars.h"
#include "TObject.h"
#include "shared.h"

#define THINNING 20

using namespace std;



//class to manage a Markov Chain Monte Carlo
//  construct with pointer to atmFitPars
class markovTools{
   public:
   
   ///////////////////
   //constructors
   markovTools(int npars);
   // This is the main constructor that is called by histoCompare
   markovTools(atmFitPars* atmpars, const char* outfilename="");
   void Init(int pars);
   void initOutputTree();

   ///////////////////////////
   //variables
   TFile* fout; //< output file
   TFile* fDiffChain;
   int nDiffSteps; //< number of events in diffChain
   TTree* diffChain;
   int nPars;  //< total number of parameters
   int nParsEffective; //< total number of non-fixed parameters;
   float effectivePars[NMCMCPARS]; //< array of non-fixed parameters
   float nominalPars[NMCMCPARS]; //< array of nominal parameters
   int  effectiveIndex[NMCMCPARS]; //< indiciies of non-fixed parameters
   int  parBin[NMCMCPARS]; //< bin of given effective parameter
   int  parIndex[NMCMCPARS]; //< atmFitPars index effective parameter
   int  parComp[NMCMCPARS]; //< MC component of given effective parameter
   int  parAtt[NMCMCPARS]; //< attribute modified by given effective parameter
   int  parIsSyst[NMCMCPARS]; //< flag = 1 if parameter is systematic (flux or xsec) parameter
   int  useDiffProposal[NMCMCPARS]; //< flag = 1 if we're using a differential MCMC to propose parameter
   int iStep;  //< counter for total step number
   float oldPars[NMCMCPARS]; //< array of parameters from previous step
   int fixPar[NMCMCPARS]; //< array of fix flags for each parameter
   float oldL; //< likelihood value of previous step;
   float tuneParameter; //< tunes the size of MCMC steps
   float varPar[NMCMCPARS]; //< stores parameter standard deviations
   float parDiff[NMCMCPARS]; //< saves the difference between mcmc steps
   int    parDiffIndex[NMCMCPARS]; //< atmFitPars index of the differential parameters
   int    ndiffpars; //< number of parameters in differential step chain
   int    nFilled; //< number of times the step tree has been filled
   int    nAccepted; //< number of times the step tree has been filled
   int    nchangethresh; //< fill threshold for changing trees
   int    noutfile;
   int    NBurnIn; //< wait this many steps before saving MCMC
   int    NThin; //< how many steps to skip before saving

   // output tree
   TTree* pathTree;
   atmFitPars* atmPars;

   /////////////////////////
   //setters
   void setFixPar(int ipar, int value){fixPar[ipar]=value;}
   void setPar(int ipar,double value){oldPars[ipar]=value;}
   void setL(double value){oldL=value; cout<<"Lset: "<<value<<endl;}
   void setParVar(int ipar,double value); //< sets parameter standard deviations
   void setDiffChain(const char* fname,int chainflg = 0); //< sets up the differential tree for DEMCMC
   void setDiffChain(TChain* ch); //< sets up the differential tree for DEMCMC
   void initDiffChain();
   ////////////////////////
   //tree setup
   void setupPathTree();

   /////////////////////////
   //MCMC Functions
   void proposeStep(double* par);  //< proposes a new step from the given parameter set
   void proposeStep(); //< propose a step from the parameters in atmFitPars
   void proposeDifferentialStep(); //< use differential chain to propose a step
   void proposePartialDiffStep(); //< use differential steps for some parameters, flat proposal for others
   int  acceptStepLnL(double newL); //< decide if step is accepted given new LnL
   int  acceptStepLnLDiff(double newL); //< decide if step is accepted given new LnL
   int acceptStep(double newL,double* par);  
   int acceptStepLnL(double newL,double* par);

   /////////////////////////
   //I/O
   void savePath();
   void setTuneParameter(double value);
   void changeFile();

   /////////////////////////
   //debugging
   void test(int itry);
   void testAtm(int itry);
   TH1D* htest;
   TH2D* htest2D;


};

#ifdef CINTMODE
#include "markovTools.cxx"
#include "globalRandom.cxx"
#endif

#endif
