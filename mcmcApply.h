#ifndef MCMCAPPLY_H
#define MCMCAPPLY_H

#include "atmFitPars.h"
#include "fqProcessedEvent.h"
#include "mcmcReader.h"

///////////////////////////////////////////////////////////////
// Class to apply parameters from MCMC for toy MC experiments
//
// This class has a pointer the attributes of some MC event (mcEvent)
// This MC event should be a PROCESSED event, processed by preProcess
//
// It also contains a pointer to an atmFitPars object (fitPars) as
// well as the MCMC cloud (mcmcPars)


class mcmcApply{
  public:

  // constructor
  mcmcApply(atmFitPars* fitpars, mcmcReader* mcmcpars, fqProcessedEvent* mcevent);
 
  // vars
  atmFitPars* fitPars;
  mcmcReader* mcmcPars;
  fqProcessedEvent* mcEvent;

  // methods
  void setFromMCMC();
  void applyPars(int iatt=-1);
  
};


#ifdef CINTMODE
#include "mcmcApply.cxx"
#endif

#endif




