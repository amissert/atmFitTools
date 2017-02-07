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
//
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

  // methods
  void setFromMCMC();
  void applyPars(int nbin, int ncomponent, float attributeTmp[], int natt);
  // reads in attributes from large array of MC points, then applies
  // the mcmc shape parameters and re-evaluates the event selection cuts
  // returns 1 for nue, 2 for numu, 0 for neither
  int applyCutsToModifiedEvent(int iev, mcLargeArray* fastevents,bool modflg=true);


};


#ifdef CINTMODE
#include "mcmcApply.cxx"
#endif

#endif




