#ifndef MCMCAPPLY_H
#define MCMCAPPLY_H

#include "atmFitPars.h"
#include "fqProcessedEvent.h"
#include "mcmcReader.h"
#include "mcLargeArray.h"
#include "eventSelectors.h"

///////////////////////////////////////////////////////////////
// Class to apply parameters from MCMC for toy MC experiments
//
// This class has a pointer the attributes of some MC event (mcEvent)
// This MC event should be a PROCESSED event (processed by preProcess)
//
// It also contains a pointer to an atmFitPars object (fitPars) as
// well as the MCMC cloud (mcmcPars)
//
// This class is useful for applying the parametesr in mcmcPars to the
// attribute[] array of a fitqun event, and then seeing if the event 
// still passes the event selection cuts defined by "eventSelectors.h"
class mcmcApply{
  public:

  // constructor
  mcmcApply(atmFitPars* fitpars, mcmcReader* mcmcpars);
 
  // vars
  atmFitPars* fitPars;
  mcmcReader* mcmcPars;
  fqProcessedEvent* mcEvent;
  int indexPIDPar;
  int indexPi0Par;
  int indexPiPPar;
  int indexRCPar;

  // best fit pars
  float bestFitPar[NMCMCPARS];


  // flags
  bool flgApplyXSecPar;
  bool flgApplyFluxPar;
  bool flgApplyNormPar;
  bool flgApplyTest;
  bool flgApplyAttSmearPar[NATTMAX];
  bool flgApplyAttBiasPar[NATTMAX];
  void setAttFlgs(int iatt, bool value){
     flgApplyAttSmearPar[iatt]=value;
     flgApplyAttBiasPar[iatt]=value;
  }
  void setAllAlphaFlags(bool value){
    flgApplyXSecPar = value;
    flgApplyFluxPar = value;
    flgApplyNormPar = value;
    return;
  }
  bool flgUseBestPar[NMCMCPARS];
  bool flgGlobalUseBestPars;


  // methods
  void setFromMCMC();
  void setAlphaFromMCMC();
  void setBetaFromMCMC();
  void setBestPars();

  // modify the reconstructed fiTQun attributes
  void applyPars(int nbin, int ncomponent, float attributeTmp[], int natt);

  // get additional weights from xsec parameters, normalization, etc
  float getXsecWeight(int mode, float enutrue);

  // get additional weights from xsec parameters, normalization, etc
  float getFluxWeight(float enutrue, int nutype);

  // test weight for debugging etc
  float getTestWeight(int nbin, int nsamp, int mode, float enutrue, int nutype=12);

  // get additional weights from xsec parameters, normalization, etc
  float getEvtWeight(int nbin, int nsamp, int mode, float enutrue, int nutype=12);
  
  // reads in attributes from large array of MC points, then applies
  // the mcmc shape parameters and re-evaluates the event selection cuts
  // returns 1 for nue, 2 for numu, 0 for neither
  int applyCutsToModifiedEvent(int iev, mcLargeArray* fastevents,bool modflg=true);

  // calculate what the best fit (mean) parameters are
  void findBestFitPars();
  void setUseBestFitSystPars(bool value);
  void setUseBestFitNormPars(bool value);


};


#ifdef CINTMODE
#include "mcmcApply.cxx"
#endif

#endif




