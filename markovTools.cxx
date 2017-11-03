#ifndef MARKOV_TOOLS
#define MARKOV_TOOLS

#include "markovTools.h"

using namespace std;



/////////////////////////////////////////////////////
// change to a new output file
void markovTools::changeFile(){

  // get current file name
  TString foutname = fout->GetName();

  // strip '.root'
  foutname.Remove(foutname.Sizeof()-6);
  // strip additional modifiers
  int nstrip = foutname.First("_");
  foutname.Remove(nstrip);

  // increment files
  noutfile++;;

  // make new file name
  TString foutnewname = foutname.Data();
  foutnewname.Append(Form("_%d",noutfile));
  foutnewname.Append(".root");

  // close old file
  cout<<"Closing file "<<fout->GetName()<<endl; 
  fout->Close();

  // make a new file
  fout = new TFile(foutnewname.Data(),"recreate");
  cout<<"Opened new file: "<<fout->GetName()<<endl;

  // setup trees, addresses, etc  
  initOutputTree(); 
 
  //
  return; 
}



/////////////////////////////////////////////////////
// set pointer to another chian
void markovTools::setDiffChain(TChain* mychain){

  diffChain = (TTree*)mychain;
  cout<<"markovTools: Using TChain with "<<diffChain->GetEntries()<<" entries!"<<endl;

  initDiffChain();

  //
  return;
}



/////////////////////////////////////////////////////
// initial setup for DEMCMC
void markovTools::initDiffChain(){

  // turn on relevant branches
  diffChain->SetBranchStatus("*",0);
  diffChain->SetBranchStatus("pardiff",1);
  diffChain->SetBranchStatus("npars",1);
  diffChain->SetBranchStatus("parindex",1);

  // set addresses
  diffChain->SetBranchAddress("pardiff",parDiff);
  diffChain->SetBranchAddress("npars",&ndiffpars);
  int atmDiffIndex[NMCMCPARS]; //< temporary array for atmFitPars indicies
  diffChain->SetBranchAddress("parindex",atmDiffIndex);

  // get number of steps
  nDiffSteps = diffChain->GetEntries();

  // set ndiffpars variable and atmDiffIndex
  diffChain->GetEntry(0);

  // by default, do not use differential proposal
  for (int ipar=0; ipar<nParsEffective; ipar++){
    useDiffProposal[ipar] = 0.;
  }

  // loop over the list of differential parameters
  for (int idiffpar = 0; idiffpar<ndiffpars; idiffpar++){
    int atmindex = atmDiffIndex[idiffpar]; //< index in full list
    // try to find a matching parameter in effective par list
    for (int imcmcpar=0; imcmcpar<nParsEffective; imcmcpar++){
      int atmindex_mcmc = parIndex[imcmcpar];
      // if a match is found, use differential chain for that parameter
      if (atmindex==atmindex_mcmc){
         useDiffProposal[imcmcpar] = 1;
         parDiffIndex[imcmcpar] = idiffpar;
         cout<<"Using differential proposal for parameter #"<<imcmcpar<<endl;
      }
    }
  }

  // ignore all but pardiff branch
  diffChain->SetBranchStatus("*",0);
  diffChain->SetBranchStatus("pardiff",1);;

  //
  return;

}


/////////////////////////////////////////////////////
// does some initial setup for differentail step tree
void markovTools::setDiffChain(const char* fname, int flgchain){

  // list of ROOT files will contain '*'
  TString names = fname;
  if (names.First("*")>=0){
    TChain* chain = new TChain("MCMCdiff");
    chain->Add(fname);
    setDiffChain(chain);
   initDiffChain();
   return;
  }

  // setup chain of par differences
  fDiffChain = new TFile(fname);
  diffChain = (TTree*)fDiffChain->Get("MCMCdiff"); 
  cout<<"markovTools: using file with tree of size "<<diffChain->GetEntries()<<endl;

  // initialize
  initDiffChain();

  // cd back to output file 
  fout->cd();
  return;
}



void markovTools::setupPathTree(){



}



void markovTools::savePath(){
  cout<<"makrovTools: path tree has "<<pathTree->GetEntries()<<" points!"<<endl;
  cout<<"markovTools: saving file "<<fout->GetName()<<endl;
  pathTree->Write();
  fout->Close();
  return;
}



void markovTools::setParVar(int ipar,double value){
  varPar[ipar] = value;
  cout<<"parameter "<<ipar<<" variance set to: "<<value<<endl;
  return;
}



void markovTools::testAtm(int ntry){

  //setup atmFitPars
  atmFitPars* fitpars = new atmFitPars(1,1,1,1);
  fitpars->setParameter(0,0.);
  fitpars->setParameter(1,0.);
  fitpars->setSysParUnc(0,1.);
  fitpars->setSysParUnc(1,1.);
  nPars = fitpars->nTotPars;
  atmPars = fitpars;

  tuneParameter = 1.0;

  int itry = 0;
  int iaccept = 0;
  double LnL;
  htest2D = new TH2D("htest2","htest2",50,-5,5,50,-5,5);
  while (itry<ntry){
    cout<<"step: "<<itry<<endl;
    double par1 = atmPars->getParameter(0);
    double par2 = atmPars->getParameter(1);
    LnL = par1*par1 + par2*par2;
    setL(LnL);
    cout<<"L: "<<oldL<<endl;
    proposeStep();
    par1 = atmPars->getParameter(0);
    par2 = atmPars->getParameter(1);
    LnL = par1*par1 + par2*par2;
    if (acceptStepLnL(LnL)){
      iaccept++;
      htest2D->Fill(atmPars->getParameter(0),atmPars->getParameter(1));
    }
    cout<<"accepted: "<<iaccept<<endl; 
    itry++;
  }
  htest2D->Draw("colz");
  //pathTree->SaveAs("pathTree.root");
  return;
}

void markovTools::test(int ntry){
  nPars = 1;
  double thePar[1] = {0.};
  int itry = 0;
  int iaccept = 0;
  htest = new TH1D("htest","htest",100,-5,5);
  while (itry<ntry){
    cout<<"step: "<<itry<<endl;
    setL(1.*(thePar[0]*thePar[0]));
    cout<<"L: "<<oldL<<endl;
    proposeStep(thePar);
    cout<<"proposal: "<<thePar[0]<<endl;
    if (acceptStepLnL(1.*(thePar[0]*thePar[0]),thePar)){
      iaccept++;
      htest->Fill(thePar[0]);
    }
    cout<<"accepted: "<<iaccept<<endl; 
    itry++;
  }
  htest->Draw();
  pathTree->SaveAs("pathTree.root");
  return;
}


/////////////////////////////////////////////
//decide if new parameters shoudl be accepted
//and write differential chain
int markovTools::acceptStepLnLDiff(double newL){

  
  // get likelihood difference
  double alpha = (oldL-newL); //< get difference in LnL
  double rand = randy->Rndm(); //< throw a random number
  int iaccept = 0; //< acceptance flag

  /////////////////////////////
  //chekc if we should accept
  if (alpha>TMath::Log(rand)){
    //accepted! new pars are now old
    for (int i=0;i<nPars;i++){
       // calculate parameter differences
       if (fixPar[i]!=1){
         parDiff[effectiveIndex[i]] = atmPars->getParameter(i) - oldPars[i];
       }
    }
    for (int i=0;i<nPars;i++){
      oldPars[i]=atmPars->getParameter(i);
      // fill array of all effective (non-fixed) parameter values
      if (fixPar[i]!=1){
        effectivePars[effectiveIndex[i]] = atmPars->getParameter(i); 
      }
      oldL = newL;
    }
#ifdef T2K
    atmPars->acceptStep();
#endif
    if ((iStep%NThin)==0) pathTree->Fill();
    iaccept = 1;
  } 
  else{
    for (int i=0;i<nPars;i++){
      //rejected, reset to old parameters
      atmPars->setParameter(i,oldPars[i]);
    }
    if ((iStep%NThin)==0) pathTree->Fill();
  }

  iStep++; //< increment global step  count
  if ((iStep%100)==0) cout<<"step: "<<iStep<<endl;

  /////////////////////////////
  return iaccept;
}


/////////////////////////////////////////////
//decide if new parameters shoudl be accepted
//this is the one used by histoCompare
int markovTools::acceptStepLnL(double newL){


  double alpha = (oldL-newL); //< get difference in LnL
  double rand = randy->Rndm(); //< throw a random number
  int iaccept = 0; //< acceptance flag
  /////////////////////////////
  //chekc if we should accept
  if (alpha>TMath::Log(rand)){
    //accepted! new pars are now old
    for (int i=0;i<nPars;i++){
      oldPars[i]=atmPars->getParameter(i);
      // fill array of all effective (non-fixed) parameter values
      if (fixPar[i]!=1){
        effectivePars[effectiveIndex[i]] = atmPars->getParameter(i); 
      }
      oldL = newL;
    }
#ifdef T2K
    atmPars->acceptStep();
#endif
    // fill output tree
    if (nAccepted>NBurnIn){
      if ((iStep%NThin)==0) pathTree->Fill();
    }
    nAccepted++;
    iaccept = 1;
  } 
  else{
    for (int i=0;i<nPars;i++){
      //rejected, reset to old parameters
      atmPars->setParameter(i,oldPars[i]);
    }
    if ((iStep%NThin)==0) pathTree->Fill();
  }
  iStep++; //< increment global step  count
  if ((iStep%100)==0) cout<<"step: "<<iStep<<endl;

  /////////////////////////////
  return iaccept;
}



int markovTools::acceptStepLnL(double newL,double* par){
  double alpha = (oldL-newL);
  cout<<"oldL: "<<oldL<<endl;
  cout<<"newL: "<<newL<<endl;
  cout<<"Ldiff: "<<alpha<<endl;
  double rand = randy->Rndm();
  cout<<"rand: "<<rand<<endl;
  int iaccept = 0;
  //cout<<"delta L: "<<alpha<<endl;
  if (alpha>TMath::Log(rand)){
  //  pathTree->SetBranchAddress("par",par);
    for (int i=0;i<nPars;i++){
      oldPars[i]=par[i];
    }
#ifdef T2K
    atmPars->acceptStep();
#endif
    if ((iStep%NThin)==0) pathTree->Fill();
    iaccept = 1;
  } 
  else{
    if ((iStep%NThin)==0) pathTree->Fill();
    for (int i=0;i<nPars;i++){
      par[i]=oldPars[i];
    }
  }
  iStep++;
  if ((iStep%100)==0) cout<<"step: "<<iStep<<endl;
  if (iaccept) cout<<" accepted! "<<endl;
  return iaccept;
}



int markovTools::acceptStep(double newL,double* par){
  double alpha = 1.;
  double lnOld;
  double lnNew;
  lnOld = TMath::Log(oldL);
  lnNew = TMath::Log(newL);
  alpha = (lnNew-lnOld);
  double rand = randy->Rndm();
  int iaccept = 0;
  if (alpha>TMath::Log(rand)){
    cout<<"accept"<<endl;
    for (int i=0;i<nPars;i++){
      oldPars[i]=par[i];
    }
    if ((iStep%NThin)==0) pathTree->Fill();
    iaccept = 1;
  } 
  else{
    cout<<"reject"<<endl;
    for (int i=0;i<nPars;i++){
      par[i]=oldPars[i];
    }
    if ((iStep%NThin)==0) pathTree->Fill();
  }
  iStep++;
  return iaccept;
}



/////////////////////////////////////////////////
//takes a poiter to a parameter array and returns
//a new array of proposed parameters
void markovTools::proposeStep(double* par){  
#ifdef T2K
  atmPars->proposeStep();
#endif
  for (int i=0;i<nPars;i++){
#ifndef T2K
    oldPars[i] = par[i];
    if (fixPar[i]!=1){
      par[i] = randy->Gaus(par[i],varPar[i]*tuneParameter);
      cout<<"Proposed parameter "<<i<<" as :"<<par[i]<<endl;
    }
#else
    if (fixPar[i]!=1) par[i] = atmPars->getPropParameter(i);
#endif
  }
  return;
}


///////////////////////////////////////////////////////////
//takes a poiter to atmFitPars and suggests new parameters
//this method uses a TTree filled with differential steps
void markovTools::proposePartialDiffStep(){  


  // get a random point in the differential chain;
  int randpoint = randy->Integer(nDiffSteps);
  diffChain->GetEntry(randpoint);

  // determins how much randomness to add
  double perturb = 0.025;

#ifndef T2K
  // loop over non-fixed pars and suggest step from chain if using it, otherwise just guess
  for (int i=0; i<nParsEffective;i++){
    // old pars are the current (stale) parameters
    int atmindex = parIndex[i];
    oldPars[atmindex] = atmPars->getParameter(atmindex);
    
    if (useDiffProposal[i]){
       int diffindex = parDiffIndex[i]; //< get index of par i in differential par list
       double epsilon = randy->Gaus(0.,varPar[atmindex])*perturb; //< random perturbation;
       epsilon = 0.;
       double diffval = parDiff[diffindex];
       double newvalue = oldPars[atmindex] + (parDiff[diffindex]*tuneParameter) + epsilon;
       atmPars->setParameter(atmindex,newvalue);
    }
    else{
      double newvalue = randy->Gaus(oldPars[atmindex],varPar[atmindex]*tuneParameter);
      atmPars->setParameter(atmindex,newvalue);
    }
  }
#else
  // need to add differential step here
  atmPars->proposeStep();
  for (int i = 0; i < nPars; ++i) {
    oldPars[i] = atmPars->getParameter(i);
  }
#endif
  return;
}



///////////////////////////////////////////////////////////
//takes a poiter to atmFitPars and suggests new parameters
//this method uses a TTree filled with differential steps
void markovTools::proposeDifferentialStep(){  

  // get a random point in the differential chain;
  int randpoint = randy->Integer(nDiffSteps);
  diffChain->GetEntry(randpoint);

  // determins how much randomness to add
  double perturb = 0.025;

#ifndef T2K
  // set atmFitPars to new values
  for (int i=0;i<nPars;i++){
    oldPars[i] = atmPars->getParameter(i);
    if (atmPars->fixPar[i]!=1){
      double epsilon = randy->Gaus(0.,varPar[i]*perturb);
      double newvalue = oldPars[i] + (parDiff[i]*tuneParameter) + perturb*randy->Gaus(0.,varPar[i]);
      atmPars->setParameter(i,newvalue);
    }
  }

#else
  // need to add differential step here
  atmPars->proposeStep();
  for (int i = 0; i < nPars; ++i) {
    oldPars[i] = atmPars->getParameter(i);
  }
#endif
  return;
}


/////////////////////////////////////////////////
//takes a poiter to atmFitPars and suggests new parameters
//
void markovTools::proposeStep(){  
#ifndef T2K
  for (int i=0;i<nPars;i++){
    oldPars[i] = atmPars->getParameter(i);
    if (atmPars->fixPar[i]!=1){
      double random = randy->Gaus(oldPars[i],atmPars->parUnc[i]*tuneParameter);
      atmPars->setParameter(i,random);
    }
  }
#else
  atmPars->proposeStep();
  for (int i = 0; i < nPars; ++i) {
    oldPars[i] = atmPars->getParameter(i);
  }
#endif
  return;
}

/////////////////////////////////////////////
// setup output tree
void markovTools::initOutputTree(){

  // make new output tree with parameters
  pathTree = new TTree("MCMCpath","MCMCpath"); //< initialize new tree for steps


  //branch setup
  pathTree->Branch("npars",&nParsEffective,"npars/I");
  pathTree->Branch("step",&iStep,"step/I");
  pathTree->Branch("par",effectivePars,"par[500]/F");
  pathTree->Branch("pardiff",parDiff,"pardiff[500]/F");
  pathTree->Branch("parnominal",nominalPars,"parnominal[500]/F");
  pathTree->Branch("parbin",parBin,"parbin[500]/I");
  pathTree->Branch("parcomp",parComp,"parcomp[500]/I");
  pathTree->Branch("paratt",parAtt,"paratt[500]/I");
  pathTree->Branch("parindex",parIndex,"parindex[500]/I");
  pathTree->Branch("parsyst",parIsSyst,"parsyst[500]/I");
  pathTree->Branch("logL",&oldL,"logL/F");

  //
  return;
}

/////////////////////////////////////////////
// initialization
void markovTools::Init(int npars){

  // output tree with parameters
  pathTree = new TTree("MCMCpath","MCMCpath"); //< initialize new tree for steps

  // set total number of parameters
  nPars = npars;

  //current step counter
  iStep=0;

  //file counter
  noutfile=0;

  //fill call counter
  nFilled=0;
  nAccepted=0;

  //counter to change files
  nchangethresh=1000000;

  //default number of burn-in steps
  NBurnIn = 0;

  //set fix arrays etc
  for (int ipar=0;ipar<nPars;ipar++){
    setParVar(ipar,atmPars->parUnc[ipar]); //< set parameter variance
    setFixPar(ipar,atmPars->fixPar[ipar]); //< set parameter fix flag
  }

  //count all non-fixed parameters
  nParsEffective = 0;
  for (int i=0; i<nPars; i++){
    // initialize parameter list with current params
    oldPars[i] = atmPars->getParameter(i);
    if (atmPars->fixPar[i]) continue;
    else{
      effectiveIndex[i] = nParsEffective;//< index used by markovTools 
      parIndex[nParsEffective] = i; //< atmFitPars index
      cout<<"parindex "<<nParsEffective<<" is "<< parIndex[nParsEffective]<<endl;
      nominalPars[nParsEffective] = atmPars->parDefaultValue[i]; //< nominal parameter value
      parBin[nParsEffective] = atmPars->getParBin(i); //< FV bin corresponding to this param
      parComp[nParsEffective] = atmPars->getParComp(i); //< MC component corresponding to this param
      parAtt[nParsEffective] = atmPars->getParAtt(i); //< fQ attribute corresponding to this param
      if (parBin[nParsEffective]>0){
        parIsSyst[nParsEffective] = 0; //< assume systematic unless changed in next block
      }
      else{
        parIsSyst[nParsEffective] = 1;
      }
      nParsEffective++;
    }
  }

  //branch setup
  pathTree->Branch("npars",&nParsEffective,"npars/I");
  pathTree->Branch("step",&iStep,"step/I");
  pathTree->Branch("par",effectivePars,"par[500]/F");
  pathTree->Branch("pardiff",parDiff,"pardiff[500]/F");
  pathTree->Branch("parnominal",nominalPars,"parnominal[500]/F");
  pathTree->Branch("parbin",parBin,"parbin[500]/I");
  pathTree->Branch("parcomp",parComp,"parcomp[500]/I");
  pathTree->Branch("paratt",parAtt,"paratt[500]/I");
  pathTree->Branch("parindex",parIndex,"parindex[500]/I");
  pathTree->Branch("parsyst",parIsSyst,"parsyst[500]/I");
  pathTree->Branch("logL",&oldL,"logL/F");

  //autosave options
  //pathTree->SetAutoSave(100);

  //done
  return;
}

/////////////////////////////////////////////
//constructor (old)
markovTools::markovTools(int npars){
  fout = new TFile("mcmctree.root","RECREATE"); //< set output file name
  pathTree = new TTree("MCMCpath","MCMCpath"); //< initialize new tree for steps
  nPars = npars; //< set total # of parameters
  iStep = 0;  //< counter of current step
  NThin = THINNING;
  
  //////////////////////////////////////////
  //branch setup
  pathTree->Branch("npars",&nPars,"npars/I");
  pathTree->Branch("step",&iStep,"step/I");
  pathTree->Branch("par",oldPars,"par[200]/D");;


}

//////////////////////////////////////////////////////////////////////////////////////
//construct from atmFitPars (this is constructor used in histoCompare so it has to work) 
markovTools::markovTools(atmFitPars* fitpars, const char* outfilename){
  TString foutname = outfilename;
  if (!foutname.CompareTo("")){
    fout = new TFile("mcmctree.root","RECREATE"); //< set output file name
  }
  else{
    fout = new TFile(foutname.Data(),"RECREATE");
  }

  cout<<"markovTools: MCMC output file: "<<fout->GetName()<<endl;

  pathTree = new TTree("MCMCpath","MCMCpath"); //< initialize new tree for steps
  atmPars = fitpars;
  
  // initialize
  Init(fitpars->nTotPars); 
  NThin = THINNING;
  // print seed
 // int randseed = randy->GetSeed();
  //cout<<"random seed: "<<randseed<<endl;

  return;
}

void markovTools::setTuneParameter(double value)
{
  tuneParameter=value;
#ifdef T2K
  atmPars->setStepSize(tuneParameter);
#endif
}



#endif
