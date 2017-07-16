#ifndef TOYMC_CXX
#define TOYMC_CXX

#include "toyMC.h"

using namespace std;




/////////////////////////////////////////////////////////////////
// apply the cuts to a modified event and see if it passes
//   returns 1 -> passed electron selection
//   returns 2 -> passed muon selection
//   returns 3 -> passed electron 1R pi selection
int toyMC::applyCutsToModifiedEvent(int iev, bool flgmod){

  // fill tmp array with "nominal" MC values
  const int natt = 4;
  float attributesTmp[natt];
  for (int iatt=0; iatt<natt; iatt++){
    attributesTmp[iatt] = fastevents->vattribute[iev][iatt];   
  }

  if (indexPIDPar>=0) cutPars.fqpid = attributesTmp[indexPIDPar];
  if (indexPi0Par>=0) cutPars.fqpi0par = attributesTmp[indexPi0Par];
  if (indexPiPPar>=0) cutPars.fqpippar = attributesTmp[indexPiPPar];
  if (indexRCPar>=0) cutPars.fqrcpar = attributesTmp[indexRCPar];

  // modify tmp array by applying the histogram shape parameters
  if (flgmod){
    // modify these cut values by the shape parameters
    modifier->applyPars(fastevents->vbin[iev],
                        fastevents->vcomponent[iev],
                        attributesTmp,
                        natt);

   // fill cut parameter structure using modified attributes
   if (indexPIDPar>=0) cutPars.fqpid = attributesTmp[indexPIDPar];
   if (indexPi0Par>=0) cutPars.fqpi0par = attributesTmp[indexPi0Par];
   if (indexPiPPar>=0) cutPars.fqpippar = attributesTmp[indexPiPPar];
   if (indexRCPar>=0) cutPars.fqrcpar = attributesTmp[indexRCPar];

  }

  // other cut pars that are not modified
  cutPars.fqmome = fastevents->vfqmumom[iev];
  cutPars.fqmommu = fastevents->vfqemom[iev];
  cutPars.nhitac = fastevents->vnhitac[iev];
  cutPars.fqnsubev = fastevents->vfqnsubev[iev];
  cutPars.fqenue = fastevents->vfqenue[iev];
  cutPars.fqenumu = fastevents->vfqenumu[iev];
  cutPars.fqnring = fastevents->vfqnring[iev];

  // see if it passes cuts
  int passnue = selectNuE(cutPars);
  int passnumu = selectNuMu(cutPars);
  int passnue1rpi = selectNuE1Rpi(cutPars);

  //
  if (passnue>0) return 1;
  if (passnumu>0) return 2;
  if (passnue1rpi>0) return 3;

  //
  return 0;
  
}



/////////////////////////////////////////////////////////////////
// Same as above but applies "core" cuts to modified event (a la
// TN 186)
//   returns:
//     -1 -> Neither core nor tail
//      0 -> Tail
//      1 -> Core -OR- e-like core for NC
//      2 -> mu-like core for NC
//      7/3: Modified for new RC parameters (always passes RC cuts)
/////////////////////////////////////////////////////////////////
int toyMC::applyCoreCutsToModifiedEvent(int iev, int nclass, bool flgmod){

  // fill tmp array with "nominal" MC values
  const int natt = 4;
  float attributesTmp[natt];
  for (int iatt=0; iatt<natt; iatt++){
    attributesTmp[iatt] = fastevents->vattribute[iev][iatt];   
  }

  // modify tmp array by applying the histogram shape parameters
  if (flgmod){
    modifier->applyPars(fastevents->vbin[iev],
                        fastevents->vcomponent[iev],
                        attributesTmp,
                        natt);
  }

  cutPars.fqpid = attributesTmp[indexPIDPar];
  cutPars.fqpi0par = attributesTmp[indexPi0Par];
  cutPars.fqpippar = attributesTmp[indexPiPPar];
  cutPars.fqrcpar = attributesTmp[indexRCPar];

  // other cut pars that are not modified
  cutPars.fqmome = fastevents->vfqmumom[iev];
  cutPars.fqmommu = fastevents->vfqemom[iev];
  cutPars.nhitac = fastevents->vnhitac[iev];
  cutPars.fqnsubev = fastevents->vfqnsubev[iev];
  cutPars.fqenue = fastevents->vfqenue[iev];
  cutPars.fqenumu = fastevents->vfqenumu[iev];
  cutPars.fqnring = fastevents->vfqnring[iev];

  // see if it passes core cuts
  // classes: 1 -> nu e CCQE
  //          2 -> nu mu CCQE
  //          3 -> nu e CCOth
  //          4 -> nu mu CCOth
  //          5 -> NC pi0
  if (nclass==1){ //< require single ring electron and CCQE
//    if (cutPars.fqpid>=0. && cutPars.fqpi0par<=0. && cutPars.fqrcpar<=0.){//< e-like, not pi0, 1R-like
    if (cutPars.fqpid>=0. && cutPars.fqpi0par<=0.){//< e-like, not pi0
      return 1;
    }
    else{
      return 0;
    }
  }


  // muon CCQE
  if (nclass==2){
//    if (cutPars.fqpid<=0 && cutPars.fqpippar<=0. && cutPars.fqrcpar<=0.){
    if (cutPars.fqpid<=0 && cutPars.fqpippar<=0.){
      return 1; // mu-lik and not pip and 1R-like so is core
    }
    else{
      return 0;
    }
  }


  // electron CCOth
  if (nclass==3){
//    if (cutPars.fqpid>=0 && cutPars.fqpi0par<=0. && cutPars.fqrcpar<=0.){
    if (cutPars.fqpid>=0 && cutPars.fqpi0par<=0.){
      return 1; // e-like and not pi0 and 1R-like
    }
    else{
      return 0;
    }
  }


  // muon CCOth
  if (nclass==4){
//    if (cutPars.fqpid<=0 && cutPars.fqpippar<=0. && cutPars.fqrcpar<=0.){
    if (cutPars.fqpid<=0 && cutPars.fqpippar<=0.){
      return 1; //< mu-lik and not pip and 1R-like
    }
    else{
      return 0.;
    }
  }
 

  // NC Pi0
  if (nclass==5){
//    if (cutPars.fqpid>=0 && cutPars.fqpi0par<=0. && cutPars.fqrcpar<=0.
//        && cutPars.fqnsubev==1 && cutPars.fqmome>100){
    if (cutPars.fqpid>=0 && cutPars.fqpi0par<=0. && cutPars.fqnsubev==1 && cutPars.fqmome>100){
      return 1; //< e-like core
    }
//    else if (cutPars.fqpid<0 && cutPars.fqpippar<=0. && cutPars.fqrcpar<=0.){
    else if (cutPars.fqpid<0 && cutPars.fqpippar<=0.){
      return 2; //< mu-like core
    }
    else{
      return 0.;
    }
  }

  // NC Pi
  if (nclass==6){
//    if (cutPars.fqpid>=0 && cutPars.fqpi0par<=0. && cutPars.fqrcpar<=0.
//        && cutPars.fqnsubev==1 && cutPars.fqmome>100){
    if (cutPars.fqpid>=0 && cutPars.fqpi0par<=0. && cutPars.fqnsubev==1 && cutPars.fqmome>100){
      return 1; //< e-like core
    }
//    else if (cutPars.fqpid<0 && cutPars.fqpippar<=0. && cutPars.fqrcpar<=0.){
    else if (cutPars.fqpid<0 && cutPars.fqpippar<=0.){
      return 2; //< mu-like core
    }
    else{
      return 0.;
    }
  }

  //
  return -1;
   
  
}



/////////////////////////////////////////////////////////////////
toyMC::toyMC(){


}



/////////////////////////////////////////////////////////////////
// set the pointers to the mc events and the post-mcmc parameters
void toyMC::setChains(TChain* chmc, TChain* chpars, int nmcevents){

  chMC = chmc;
  chPars = chpars;

  mcEvent = new fqProcessedEvent(chMC);
  fastevents = new mcLargeArray(chMC,nmcevents);
  nMCevents = nmcevents;

  mcmcPars = new mcmcReader(chPars);

  return;
}



////////////////////////////////////////////////////////////////
//  read parameters from a random point in the MCMC cloud
int toyMC::getRandomMCMCPoint(){
  int nmax = chPars->GetEntries();
  int randpoint = randy->Integer(nmax);
  cout<<"getting mcmc point: "<<randpoint<<endl;
  chPars->GetEntry(randpoint);
  return randpoint;
}



////////////////////////////////////////////////////////////////
// get thec combined uncertainties for all events
void toyMC::makeCombinedUncertainty(int nmcmcpts){

  // setup containter for t2k sample
  t2kToys = new t2kSample("_toymc",1,1);

  // get list of mc events
  int nevmax = chMC->GetEntries();
  if (nMCevents>nevmax) nMCevents = nevmax;

  // get list of points in mcmc parameter space
  cout<<"Making list of MCMC points"<<endl;
  randomList* mcmclist = new randomList(nmcmcpts,chPars->GetEntries(),nmcmcpts);

  // loop over mcmc points
  for (int i=0; i<nmcmcpts; i++){

    // read in shape parameters
    cout<<"getting event"<<mcmclist->getAt(i)<<endl;
    chPars->GetEntry(mcmclist->getAt(i));

    // modify attributes using thes parameters
    modifier->setFromMCMC();

    // loop over T2K MC events
    for (int iev=0; iev<nMCevents; iev++){
     
      // apply cuts
      int ipass = applyCutsToModifiedEvent(iev);
      if (ipass==0) continue;

      // fill histos
      if (ipass==1) t2kToys->hEnuElectron->Fill(fastevents->vfqenue[iev], fastevents->vweight[iev]);
      if (ipass==2) t2kToys->hEnuMuon->Fill(fastevents->vfqenumu[iev], fastevents->vweight[iev]);
    }
    t2kToys->finishToyRun();
  }
  
  t2kToys->calcUncertainties();

  return;

}



////////////////////////////////////////////////////////////////
// fill the SKError class
void toyMC::fillMarginalizedSKErr(const int ntoys, const int nmarg, int nbinning, int flgcustom, int effdef){

  // make error container
  skErr = new SKError(ntoys);
  skErr->initHistos(nbinning);
  skErr->effDefinition = effdef;
  skErr->Nmarginal = nmarg;

  // decide which parameters to use
  if (!flgcustom){
    modifier->flgApplyXSecPar = true;
    modifier->flgApplyFluxPar = true;
    modifier->flgApplyNormPar= true;
  }

  // start with parameters at default values
  modifier->fitPars->resetDefaults();

  // reset histogram contents
  skErr->resetHistos();

  // setup mcmc samples
  int ntoys_total = ntoys+nmarg;
  uniformList* mcmc_list = new uniformList(ntoys_total,0,chPars->GetEntries());
  int mcmc_list_beta[ntoys]; 
  int mcmc_list_alpha[nmarg]; 
  for ( int i=0 ; i<ntoys; i++ ) {
    mcmc_list_beta[i] = mcmc_list->getAt(i);
  }
  for ( int i=ntoys ; i<ntoys+nmarg; i++ ) {
    mcmc_list_alpha[i] = mcmc_list->getAt(i);
  }

  // determine max MC events
  int nevmax = chMC->GetEntries();
  if (nMCevents>nevmax) nMCevents = nevmax;


  // fill nominal values///////////////////////////////////////////////////
  // loop over T2K MC events
  for (int iev=0; iev<nMCevents; iev++){

    // get true MC event class for event "iev"
    int nclass = skErr->getClassMC(fastevents->vnutype[iev],
                                     fastevents->vmode[iev],
                                     fastevents->vcomponent[iev],
                                     fastevents->vfqemom[iev],
                                     fastevents->vfqnsubev[iev],
                                     fastevents->vfqtowall[iev],
                                     fastevents->vfqwall[iev]);
    
    // passes core selection?
    int iscore = applyCoreCutsToModifiedEvent(iev,nclass,true);

    // if its core or tail fill histos
    if (iscore<0) continue; //< skip unclassified events

    // get the new event weight
    float   ww = fastevents->vweight[iev]
                 *modifier->getEvtWeight( fastevents->vbin[iev], fastevents->vsample[iev]
                                         ,fastevents->vmode[iev], fastevents->vpmomv[iev]
                                         ,fastevents->vnutype[iev]);
     

    // fill total histos
    // the "true" flag means we add it to the "total" evis histogram for this class
    skErr->addEventBase(nclass,fastevents->vfqemom[iev],ww);

    if (iscore>0){
      // fill core histos
    // the "true" flag means we add it to the "core" evis histogram for this class
      skErr->addEventCore(nclass,fastevents->vfqemom[iev],ww,iscore);
    }
     
  }
  skErr->addNominal();



  // fill marginalized SK errors values///////////////////////////////////////////////////
  for (int j=0; j<ntoys; j++){

    // set beta parameters
    cout<<"working on toy: "<<j<<endl;
    chPars->GetEntry(mcmc_list_beta[j]);
    modifier->setBetaFromMCMC();
    modifier->setBestPars();

    for (int i=0; i<nmarg; i++){

      // set alpha parameters
      chPars->GetEntry(mcmc_list_alpha[i]);
      modifier->setAlphaFromMCMC();
      modifier->setBestPars();

      // reset histogram contents
      skErr->resetHistos();

      // loop over T2K MC events
      for (int iev=0; iev<nMCevents; iev++){


        // get true MC event class for event "iev"
       // int nclass = skErr->getClassMC(fastevents->vnutype[iev],
       //                                fastevents->vmode[iev],
       //                                fastevents->vcomponent[iev],
       //                                fastevents->vfqemom[iev],
       //                                fastevents->vfqnsubev[iev],
       //                                fastevents->vfqtowall[iev],
       //                                fastevents->vfqwall[iev]);
       int nclass = skErr->getClassMC(fastevents->vnutype[iev],
                                       fastevents->vmode[iev],
                                       fastevents->vcomponent[iev],
                                       fastevents->vfqemom[iev],
                                       fastevents->vfqnsubev[iev],
                                       fastevents->vfqtowall[iev],
                                       fastevents->vfqwall[iev]);
    
        // passes core selection?
        int iscore = applyCoreCutsToModifiedEvent(iev,nclass,true);

        // if its core or tail fill histos
        if (iscore < 0) continue; //< skip unclassified events

        // get the new event weight
        float   ww = fastevents->vweight[iev]
                   *modifier->getEvtWeight( fastevents->vbin[iev], fastevents->vsample[iev]
                                           ,fastevents->vmode[iev], fastevents->vpmomv[iev]
                                           ,fastevents->vnutype[iev]); 

        // fill total histos
//        skErr->addEvent(nclass,fastevents->vfqemom[iev],ww,true);
        skErr->addEventBase(nclass,fastevents->vfqemom[iev],ww);
        if (iscore>0){
          // fill core histos
//          skErr->addEvent(nclass,fastevents->vfqemom[iev],ww,false);
          skErr->addEventCore(nclass,fastevents->vfqemom[iev],ww,iscore);
        }

      } //< MC eventloop

      // count the events in each class
      skErr->fillNeventArrays(i);

    } //< marginalization loop 

    // save marginalized efficiency
    skErr->marginalize(j);

  }

  skErr->calcCov();
  skErr->calcErrors();
  skErr->printErrors();

  //
  return;
}




////////////////////////////////////////////////////////////////
// fill the SKError class
void toyMC::fillSKErrors(int ntoys,int nbinning, int flgcustom, int effdef){

  // ring-counting map
  RCMap* rcMap = new RCMap("./data/RCMap.root");

  // make error container
  skErr = new SKError(ntoys);
  skErr->initHistos(nbinning);

  // decide which parameters to use
  if (!flgcustom){
    modifier->flgApplyXSecPar = true;
    modifier->flgApplyFluxPar = true;
    modifier->flgApplyNormPar= true;
  }

  // start with parameters at default values
  modifier->fitPars->resetDefaults();

  // get list of random MCMC points from the MCMC point file
  int nmcmcmax = chPars->GetEntries();
  uniformList* mcmc_points = new uniformList(ntoys,0,chPars->GetEntries());

  // determine max events
  int nevmax = chMC->GetEntries();
  if (nMCevents>nevmax) nMCevents = nevmax;

  // reset histogram contents
  skErr->resetHistos();

  // fill nominal values
  for (int iev=0; iev<nMCevents; iev++){

    // get true MC event class for event "iev"
    int nclass = skErr->getClassMC(fastevents->vnutype[iev],
                                     fastevents->vmode[iev],
                                     fastevents->vcomponent[iev],
                                     fastevents->vfqemom[iev],
                                     fastevents->vfqnsubev[iev],
                                     fastevents->vfqtowall[iev],
                                     fastevents->vfqwall[iev]);
    if (nclass==0) continue;

    // passes core selection?
    int iscore = applyCoreCutsToModifiedEvent(iev,nclass,false);

    // if it's core or tail fill histos
    if (iscore<0) continue; //< skip unclassified events

    // get the new event weight
    float   ww = fastevents->vweight[iev]
                 *modifier->getEvtWeight( fastevents->vbin[iev], fastevents->vsample[iev]
                                         ,fastevents->vmode[iev], fastevents->vpmomv[iev]
                                         ,fastevents->vnutype[iev]);
    if (flgRCWeight){
      ww *= rcMap->get1RWeight(cutPars.fqrcpar);
    }

    // the "true" flag means we add it to the "total" evis histogram for this class
    skErr->addEventBase(nclass,fastevents->vfqemom[iev],ww);
    if (iscore>0){
      skErr->addEventCore(nclass,fastevents->vfqemom[iev],ww,iscore);
    }
     
  }
  skErr->addNominal();


  // toyMC:loop over mcmc points
  int itoy=0;
  while (true){

    // read in fit parameters
    
    int ientry = mcmc_points->getAt(itoy);
    cout<<"getting event"<<ientry<<endl;
    if (ientry>=nmcmcmax) break; 
    chPars->GetEntry(ientry);

    // modify attributes using thes parameters
    modifier->setFromMCMC();
    modifier->setBestPars();

    // reset histogram contents
    skErr->resetHistos();

    // loop over T2K MC events
    for (int iev=0; iev<nMCevents; iev++){

      // get true MC event class for event "iev"
      int nclass = skErr->getClassMC(fastevents->vnutype[iev],
                                     fastevents->vmode[iev],
                                     fastevents->vcomponent[iev],
                                     fastevents->vfqemom[iev],
                                     fastevents->vfqnsubev[iev],
                                     fastevents->vfqtowall[iev],
                                     fastevents->vfqwall[iev]);
    
      if (nclass==0) continue;

      // passes core selection?
      int iscore = applyCoreCutsToModifiedEvent(iev,nclass,true);

      // if its core or tail fill histos
      if (iscore<0) continue; //< skip unclassified events

      // get the new event weight
      float   ww = fastevents->vweight[iev]
                   *modifier->getEvtWeight( fastevents->vbin[iev], fastevents->vsample[iev]
                                           ,fastevents->vmode[iev], fastevents->vpmomv[iev]
                                           ,fastevents->vnutype[iev]);
      if (flgRCWeight){
        ww *= rcMap->get1RWeight(cutPars.fqrcpar);
      }

      // fill total histos
      skErr->addEventBase(nclass,fastevents->vfqemom[iev],ww);
      if (iscore>0){
        skErr->addEventCore(nclass,fastevents->vfqemom[iev],ww,iscore);
      }
     
    }

    // save toy histo contents
    skErr->addToy(itoy);
    itoy++;
  }
  skErr->Ntoys = itoy;

  skErr->calcCov();
  skErr->calcErrors();
  skErr->printErrors();

  cout<<"filled "<<itoy<<" toys!"<<endl;
  //
  return;
}



////////////////////////////////////////////////////////////////
// Get event catagory. Current catagories:
//   1 -> CCQE
//   2 -> CCnQE
//   3 -> CCMisID
//   4 -> NC
//   0 -> Uncatagorized
int toyMC::getEventCatagory(int iev, int inutype){
      

      // is NC?
      if (fastevents->vmode[iev]>=30){
        return 4;
      }

      //  CC?
      
      // mid-IDed
      if (fastevents->vnutype[iev]!=inutype){
        return 3;
      }
      
      // ccqe
      else if (fastevents->vmode[iev]==1) {
        return 1;
      }

      else {
        return 2;
      }
      
      //
      return 0;
}



////////////////////////////////////////////////////////////////
// make map of uncertainties in different (wall,towall) regions
void toyMC::makeFVUncMap(int nmcmcpts, int nselection, const char* outfile, int fvbintype){

  // get name of selection
  TString selection_name;
  int selection_nutype;
  if (nselection==1){
    selection_name = Form("NuE_Bintype%d",fvbintype);
    selection_nutype = 12;
  }
  else if (nselection==3){
    selection_name = Form("NuE_1Pi_Bintype%d",fvbintype);
    selection_nutype = 12;
  }
  else if (nselection==2){
    selection_name = Form("NuMu_Bintype%d",fvbintype);
    selection_nutype = 14;
  }

  // make array of histos
  cout<<"makeFVUncMap: Initializing array of histograms..."<<endl;
  TH1D* hE; //< nu energy binning...this is passed on optimusPrime
  if (nselection==2){
     hE = new TH1D(selection_name.Data(),selection_name.Data(),EnuNBins,EnuBinning);
  }
  else if (nselection==1 || nselection==3){
     hE = new TH1D(selection_name.Data(),selection_name.Data(),EnuNBinsElectron,EnuBinningElectron);
  }

  // make FV histogram for finding bins
  // can use different bin type flags
  cout<<"makeFVUncMap: Initializing seed for FV histogram arrays..."<<endl;
  TH2FV* hfv = new TH2FV(Form("hfv_bintype%d",fvbintype),fvbintype);
  hfv->Draw();

  // array of nu energy and FV histograms to be filled
  cout<<"makeFVUncMap: Initializing histogram arrays"<<endl;
  modHistoArrayFV* hArrFV = new modHistoArrayFV(hE,hfv,nmcmcpts);

  // get list of mc events
  int nevmax = chMC->GetEntries();
  if (nMCevents>nevmax) nMCevents = nevmax;

  // get list of points in mcmc parameter space
  cout<<"Making list of MCMC points"<<endl;
  randomList* mcmclist = new randomList(nmcmcpts,chPars->GetEntries(),nmcmcpts);

  // start with parameters at default values
  modifier->fitPars->resetDefaults();

  // loop over mcmc points
  for (int i=0; i<nmcmcpts; i++){

    // read in parameters
    cout<<"getting event"<<mcmclist->getAt(i)<<endl;
    chPars->GetEntry(mcmclist->getAt(i));

    // modify attributes using thes parameters
    modifier->setFromMCMC();

    // loop over T2K MC events
    for (int iev=0; iev<nMCevents; iev++){

      // apply parameters and see if it passes cuts
      int ipass = 0;
      float ww = 1.0;
      if (i!=0){
        ipass = applyCutsToModifiedEvent(iev,true); 
        ww *= modifier->getEvtWeight( fastevents->vbin[iev], fastevents->vsample[iev], fastevents->vmode[iev], fastevents->vpmomv[iev] );
      }
      else{
        ipass = applyCutsToModifiedEvent(iev,false);
      }

      // if it passes fill histos
      if (ipass!=nselection) continue;

      // modified nu energy
      float enu = fastevents->vfqenumu[iev];

      // fill total nev
      int fvbin = hArrFV->hFV[i]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],ww*fastevents->vweight[iev]) - 1;

      // fill FV histos for different event catagories
      //   1 -> CCQE
      //   2 -> CCnQE
      //   3 -> CCMisID
      //   4 -> NC
      //   0 -> Uncatagorized
      int catagory = getEventCatagory(iev,selection_nutype);
      // fill histograms based on catagory
      if (catagory==1){ hArrFV->hFVCCQE[i]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],ww*fastevents->vweight[iev]);}
      if (catagory==2){ hArrFV->hFVCCnQE[i]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],ww*fastevents->vweight[iev]);}
      if (catagory==3){ hArrFV->hFVCCWrong[i]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],ww*fastevents->vweight[iev]);}
      if (catagory==4){ hArrFV->hFVNC[i]->Fill(fastevents->vfqtowall[iev],fastevents->vfqwall[iev],ww*fastevents->vweight[iev]);}
      if (catagory==0){cout<<"event # "<<iev<<" has not a catagory..."<<endl;}
      if (fvbin>=0) hArrFV->getHistogram(i,fvbin)->Fill(enu, ww*fastevents->vweight[iev]);

    }
  }

  // calculate summary and save output
  hArrFV->calcSummary();
  hArrFV->saveSummary(outfile);
  hArrFV->saveClose();
  return;



}




void toyMC::setAtmFitPars(const char* parfile){
  
  fitPars = new atmFitPars(parfile); 

  modifier = new mcmcApply(fitPars, mcmcPars);

}



////////////////////////////////////////////////////////////////
void toyMC::setCompare(histoCompare* hc){

 hCompare = hc;

 modifier = new mcmcApply(hCompare->thePars, mcmcPars);

 return;
}

















#endif
