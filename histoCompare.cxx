#ifndef HISTOCOMPARE_C
#define HISTOCOMPARE_C

#include "histoCompare.h"


histoCompare* histoCompare::staticthis;

//////////////////////////////////////////////////
// Print summary histograms after running the fit
void histoCompare::printFitSummary(const char* outdir){

  // setup new canvas
  TCanvas *can = new TCanvas("can","can",700,900);
  can->Divide(2,3);


  // plot root name
  TString plotroot = "fit_result_";

  // loop over attributes, bins, and samples and print all 
  for (int isamp=0; isamp<nSamp; isamp++){
    for (int iatt=0; iatt<nAtt; iatt++){

       // get plot name
       TString plotname = outdir;
       plotname.Append(plotroot.Data());
       plotname.Append(Form("sample_%d_",isamp));
       plotname.Append(Form("attribute_%d",iatt));
       plotname.Append(".png");

       // fill each pad with the fit results
       for (int ibin=0; ibin<nBin; ibin++){
         can->cd(ibin+1);
         showFitResult(isamp,ibin,iatt);
       }
       
       cout<<"Saving plot: "<<plotname.Data()<<endl;
       can->Print(plotname.Data());
    }
  }

  return;
}



/////////////////////////////////////////////
//get a rough estimation (fron 1D likelihood proflie)
//of the uncertainty of each fit parameter
//Requires a guess of parUnc[] in atmFitPars initialization
//This should be run before running MCMC, since these
//rough uncertainties are used to proposes MCMC steps
void histoCompare::calcRoughParErr(){
   cout<<"histoCompare: "<<"calculating rough uncertainties"<<endl;  
  //print final results
  for (int ipar=0;ipar<thePars->nTotPars;ipar++){
    if (thePars->fixPar[ipar]) continue;
    errParLo[ipar]=getErrLo(ipar);
    errParHi[ipar]=getErrHi(ipar);
    double errnew = errParHi[ipar]-errParLo[ipar];
    if (errnew>0.) thePars->parUnc[ipar]=(errParHi[ipar]-errParLo[ipar]);
    cout<<"  PAR "<<ipar<<" FIT RESULT: "<<thePars->pars[ipar]<<" +/- : "<<thePars->parUnc[ipar]<<endl;
  }
  return;
}

void histoCompare::readFitPars(const char* filename){
  thePars->readPars(filename);
  return;
}

void histoCompare::timetest(int ntry){
  clock_t t1,t2;
  double par[100];
  int parindex = 0;
  for (int ibin=0;ibin<nBin;ibin++){
    for (int icomp=0;icomp<nComp;icomp++){
      for (int iatt=0;iatt<nAtt;iatt++){
        for (int imod=0;imod<2;imod++){
           par[parindex]=Par[ibin][icomp][iatt][imod];
           parindex++;
        }
      }
    }
  }
  double diff;
  int itry=0;
  t1=clock();
  while (itry<ntry){
    //getTotLnL1D(result,par); 
    getTotLnL();
    itry++;
  }
  t2=clock();
  diff = ((double)t2-(double)t1)/(double)ntry;
  double diff1 = diff;
  cout<<"time 1: "<<diff<<endl;
  itry=0;
  t1=clock();
  while (itry<ntry){
    //getTotLnL1D(result,par); 
    getTotLnL();
    itry++;
  }
  t2=clock();
  diff = ((double)t2-(double)t1)/(double)ntry;
  cout<<"time 2: "<<diff<<endl;
  cout<<"changes: "<<diff1-diff<<endl;
  return;
}


void histoCompare::tuneDEMCMC(int ncycles,int nsteps, double goal){
 
  double result = 0.;
//  markovTools* mcmc = new markovTools(thePars); //< create markovTools object
  mcmctools = new markovTools(thePars); //< create markovTools object
  markovTools* mcmc = mcmctools;
  mcmc->setTuneParameter(tunePar);
  mcmc->setDiffChain(diffChainFileName.Data());

  //set initial state
  result = getTotLnL();
  mcmc->setL(result);//< sets the initial likelihood
  double Linit = result;

  //run tuning
  for (int icycle=0;icycle<ncycles;icycle++){
    double xaccepted=0.;
    int   istep=0;
    while (istep<nsteps){
      cout<<"----------"<<"step : "<<istep<<"---------------"<<endl;;
//      mcmc->proposeDifferentialStep(); //< propose a new step using differential mcmc
      mcmc->proposePartialDiffStep(); //< propose a new step using differential mcmc
      result = getTotLnL();
      cout<<"hc: Likelihood "<<Linit<<" -> "<<result<<" diff: "<<result-Linit<<endl;
      if (mcmc->acceptStepLnL(result)){ //< check if new step is accepted
        xaccepted++; 
      }
      istep = mcmc->iStep;
    }

    double rate = xaccepted/(double)nsteps;
    cout<<"acceptance rate: "<<rate<<endl;
    cout<<"tune parameter: "<<tunePar<<endl;
    if ((rate>22.0)&&(rate<30.0)) return;
    tunePar*=(rate/goal);
    cout<<"new tune parameter: "<<tunePar<<endl;;
  }
  return; 
}



void histoCompare::tuneMCMC(int ncycles,int nsteps,double goal){
 
  double result = 0.;
  markovTools* mcmc = new markovTools(thePars); //< create markovTools object
  mcmc->setTuneParameter(tunePar);


  //set initial state
  result = getTotLnL();
  mcmc->setL(result);//< sets the initial likelihood
  double Linit = result;

  //run tuning
  for (int icycle=0;icycle<ncycles;icycle++){
    double xaccepted=0.;
    int   istep=0;
    double Lprevstep = Linit;
    double Lprevaccept = Linit;
    while (istep<nsteps){
      cout<<"----------"<<"step : "<<istep<<"---------------"<<endl;;
      mcmc->proposeStep(); //< propose a new step
      result = getTotLnL();
      cout<<"hc: Likelihood "<<Lprevstep<<" -> "<<result<<" diff: "<<result-Lprevstep<<endl;
      cout<<"hc: From last step: "<<result - Lprevaccept<<endl;
      if (mcmc->acceptStepLnL(result)){ //< check if new step is accepted
        cout<<"ACCEPTED!"<<" total distance: "<<result-Linit<<endl;
        Lprevaccept = result;
        xaccepted++; 
      }
      istep = mcmc->iStep;
      Lprevstep = result;
    }

    double rate = xaccepted/(double)nsteps;
    cout<<"acceptance rate: "<<rate<<endl;
    cout<<"tune parameter: "<<tunePar<<endl;
    tunePar*=(rate/goal);
    cout<<"new tune parameter: "<<tunePar<<endl;;
  }
  return; 
}







///////////////////////////////////////////////
//Run a MCMC of length nsteps
void histoCompare::runDEMCMC(int nsteps){

  ///////////////////////////////////////
  //setup mcmc tools
  markovTools* mc = new markovTools(thePars,MCMCOutputFile.Data());
  mc->setTuneParameter(tunePar);

  //see if we're using a chain or a single file
  if (diffChainFileName.First("*")>=0){
    TChain* mychain = new TChain("MCMCdiff");
    mychain->Add(diffChainFileName.Data());
    mc->setDiffChain(mychain);
  }
  else{
   mc->setDiffChain(diffChainFileName.Data());
  }

  ///////////////////////////////////////////////
  //if no steps spefied, use number from par file
  if (nsteps<0){
    nsteps = MCMCNSteps;
  }

  ///////////////////////////////////////////////
  //set burn-in before parameters are saved
  mc->NBurnIn = MCMCNBurnIn;

  ///////////////////////////////////////////////
  //fill parameter array and set uncertainties
  for (int ipar=0;ipar<thePars->nTotPars;ipar++){
    mc->setParVar(ipar,thePars->parUnc[ipar]);
  }
  
  ///////////////////////////////////////////////////
  //set initial state
  double result = getTotLnL();
  mc->setL(result);//< sets the initial likelihood

  //loop through steps
  int currentstep=0;;
  while (currentstep<nsteps){
    currentstep = mc->iStep;
    mc->proposePartialDiffStep();
    result = getTotLnL();
    mc->acceptStepLnL(result);
  }

  ///////////////////////
  //save results
  mc->savePath();

  //////////////////////////
  return;
}



///////////////////////////////////////////////
//Run a MCMC of length nsteps
void histoCompare::runMCMC(int nsteps){


  ///////////////////////////////////////////////
  // print seed
  int randseed = randy->GetSeed();
  cout<<"random seed: "<<randseed<<endl;

  ///////////////////////////////////////////////
  //if no steps spefied, use number from par file
  if (nsteps<0){
    nsteps = MCMCNSteps;
  }

  ///////////////////////////////////////
  //setup mcmc tools
  markovTools* mc = new markovTools(thePars,MCMCOutputFile.Data());
  mc->setTuneParameter(tunePar);

  ///////////////////////////////////////////////
  //fill parameter array and set uncertainties
  for (int ipar=0;ipar<thePars->nTotPars;ipar++){
    mc->setParVar(ipar,thePars->parUnc[ipar]);
  }
  
  ///////////////////////////////////////////////////
  //set initial state
  double result = getTotLnL();
  mc->setL(result);//< sets the initial likelihood

  ///////////////////////////////////////////////
  //Number of burn-in steps
  mc->NBurnIn = MCMCNBurnIn;

  ///////////////////////////////////////////////
  // print tune parameter
  cout<<"MCMC tune parameter: "<<tunePar<<endl;


  //loop through steps
  int currentstep=0;
  while (currentstep<nsteps){
    currentstep = mc->iStep;
    //mc->proposeStep(par); //< fills par[] with proposed params
    mc->proposeStep();
    //getTotLnL1D(result, npars,par);   
    result = getTotLnL();
    //mc->acceptStepLnL(result,par); //< if step is accepted, istep++, and written to histo
    int isaccept = mc->acceptStepLnL(result);
  }

  ///////////////////////
  //save results
  mc->savePath();

  //////////////////////////
  return;
}



void histoCompare::saveFitPars(const char* filename){
  thePars->savePars(filename);
  return;
}



double histoCompare::getErrHi(int ipar){

  //scan through log likelihood by decreasing parameter until threshold is reached
  double Lthresh = 1.0; //likelihood threshold
  double precision = 0.2;
  double Ldiff = 0.;  //difference in likelihood from current value
  double Lbest = getTotLnL(); //current likelihood value
  double parbest = thePars->pars[ipar];
  double parval = thePars->pars[ipar]; 
  double dpar = thePars->parUnc[ipar];
  double hierr;
  int ntry = 0;
  int ntrymax;


  //guess at a good dpar
  
  thePars->setParameter(ipar,parval+(3.*dpar));
  Ldiff += TMath::Abs(Lbest-getTotLnL());
  if (Ldiff>0.){
    dpar = (Lthresh*3*dpar)/(4*Ldiff);
  }
  thePars->setParameter(ipar,parbest);

  thePars->setParameter(ipar,parval+(4.*dpar));
  Ldiff = TMath::Abs(Lbest-getTotLnL());
  thePars->setParameter(ipar,parval-(4.*dpar));
  Ldiff += TMath::Abs(Lbest-getTotLnL());
  if (Ldiff>0.){
    dpar = (Lthresh*dpar*8)/(8.*Ldiff);
  }
  thePars->setParameter(ipar,parbest);

  //coarse search
  ntrymax = 10;
  Ldiff = 0.;
  while (Ldiff<Lthresh){
    parval+=dpar;
    thePars->setParameter(ipar,parval); //modify parameter
    Ldiff = TMath::Abs(Lbest-getTotLnL()); //check L difference;
    ntry++;
    if (ntry>ntrymax){
      break;
    }
  }
  parval-=dpar; //< go back one step
  thePars->setParameter(ipar,parval);

  //fine search
  Ldiff = 0;
  dpar*=0.10;
  ntry=0;
  ntrymax = 50;

  while ((Ldiff<1)){
    parval+=dpar;
    thePars->setParameter(ipar,parval); //modify paramete
    Ldiff = fabs(Lbest-getTotLnL()); //check L difference
    ntry++;
    if (ntry>ntrymax) break;
  }
  parval-=dpar;
  thePars->setParameter(ipar,parval);
  hierr = thePars->pars[ipar];
  thePars->setParameter(ipar,parbest);
  return hierr; 
}




double histoCompare::getErrLo(int ipar){

  //scan through log likelihood by decreasing parameter until threshold is reached
  double Lthresh = 1.0; //likelihood threshold
  double precision = 0.2;
  double Ldiff = 0.;  //difference in likelihood from current value
  double Lbest = getTotLnL(); //current likelihood value
  double parbest = thePars->pars[ipar];
  double parval = thePars->pars[ipar]; 
  double dpar = thePars->parUnc[ipar];
  double loerr;
  int ntry = 0;
  int ntrymax;


  //guess at a good dpar
  thePars->setParameter(ipar,parval-(3.*dpar));
  Ldiff += TMath::Abs(Lbest-getTotLnL());
  if (Ldiff>0.){
    dpar = (Lthresh*3*dpar)/(4*Ldiff);
  }
  thePars->setParameter(ipar,parbest);

  //coarse search
  ntrymax = 10;
  while (Ldiff<Lthresh){
    parval-=dpar;
    thePars->setParameter(ipar,parval); //modify parameter
    Ldiff = TMath::Abs(Lbest-getTotLnL()); //check L difference;
    ntry++;
    if (ntry>ntrymax){
      break;
    }
  }
  parval+=dpar; //< go back one step
  thePars->setParameter(ipar,parval);

  //fine search
  Ldiff = 0;
  dpar*=0.10;
  ntry=0;
  ntrymax = 50;
  while ((Ldiff<1)){
    parval-=dpar;
    thePars->setParameter(ipar,parval); //modify paramete
    Ldiff = fabs(Lbest-getTotLnL()); //check L difference
    ntry++;
    if (ntry>ntrymax) break;
  }
  parval+=dpar;
  thePars->setParameter(ipar,parval);
  loerr = thePars->pars[ipar];
  thePars->setParameter(ipar,parbest);
  return loerr; 
}



void histoCompare::profileL(int ipar, double range, int npts, int sameflg){
  TString pname = "p";
  pname.Append("_profile.png");
  int iprofile =0;
  if (sameflg) iprofile = 1;
  double bestpoint = thePars->pars[ipar];
  cout<<"Fitted value: "<<bestpoint<<endl;
  double dx = range/(double)npts;
  double xx = bestpoint - (range/2.);
  double ll;
  double lbest = getTotLnL();
  if (hProf[iprofile]!=NULL) hProf[iprofile]->Delete();
  hProf[iprofile] = new TH1D("hprof","hprof",npts,xx,(xx+range));
  cout<<"histoCompare: profiling parameter "<<ipar<<endl;
  for (int ipoint=0;ipoint<npts;ipoint++){
    thePars->setParameter(ipar,xx);
    ll = getTotLnL();
    hProf[iprofile]->SetBinContent(ipoint+1,ll-lbest);
    xx+=dx;
  } 
  //set parameter back to initial value
  thePars->setParameter(ipar,bestpoint);
  hProf[iprofile]->SetLineColor(9);
  hProf[iprofile]->GetYaxis()->SetTitle("#Delta #chi^{2}");
  TString xname = "parameter ";
  xname.Append(Form("%d",ipar));
  hProf[iprofile]->GetXaxis()->SetTitle(xname.Data());
  hProf[iprofile]->Draw("c");
  if (sameflg) hProf[iprofile]->SetLineColor(kRed);
  if (sameflg){ 
    hProf[0]->Draw("c");
    hProf[iprofile]->Draw("samec");
  }
  else{
    hProf[0]->Draw("c");
  }
  cc->Print(pname.Data());
  return;
}



void histoCompare::profileL(int ibin, int icomp, int iatt, int imod, double range, int npts){
  TString pname = "p";
  pname.Append("_profile.png");
  double bestpoint = thePars->getHistoParameter(ibin,icomp,iatt,imod);


  cout<<"Fitted value: "<<bestpoint<<endl;
  double dx = range/(double)npts;
  double xx = bestpoint - (range/2.);
  double ll;
  double lbest = getTotLnL();
  if (hProf[0]) hProf[0]->Delete();
  hProf[0] = new TH1D("hprof","hprof",npts,xx,(xx+range));
  for (int ipoint=0;ipoint<npts;ipoint++){
    thePars->setParameter(ibin,icomp,iatt,imod,xx);
    ll = getTotLnL();
    hProf[0]->SetBinContent(ipoint+1,ll-lbest);

    xx+=dx;
  } 
  //set parameter back to initial value
  thePars->setParameter(ibin,icomp,iatt,imod,bestpoint);
  hProf[0]->SetLineColor(9);
  hProf[0]->Draw("c");
  cc->Print(pname.Data());
  return;
}

void histoCompare::showFitPars(int ibin,int iatt,int imod){
  const int nbinstot = nComp;
  double X[nbinstot];
  double EXL[nbinstot];
  double EXH[nbinstot];
  double EYL[nbinstot];
  double EYH[nbinstot];
  double Y[nbinstot];
  int parindex;
  if (gPar) gPar->Delete();
  hPar  = new TH1D("hpar","hpar",nbinstot,0,nbinstot);
  hParErrLo = new TH1D("hparerrlo","hparerrlo",nbinstot,0,nbinstot);
  hParErrHi = new TH1D("hparerrhi","hparerrhi",nbinstot,0,nbinstot);
  for (int icomp=0;icomp<nComp;icomp++){
    parindex=thePars->getParIndex(ibin,icomp,iatt,imod);
    hPar->SetBinContent(icomp+1,thePars->getParameter(parindex));
    X[icomp]=(double)hPar->GetBinCenter(icomp+1);
    Y[icomp]=(double)thePars->getParameter(parindex);
    EYL[icomp]=(double)thePars->getParameter(parindex) -errParLo[parindex];
    EYH[icomp]=errParHi[parindex] - (double)thePars->getParameter(parindex);
    hParErrLo->SetBinContent(icomp+1,(double)EYL[icomp]);
    hParErrHi->SetBinContent(icomp+1,(double)EYH[icomp]);
    EXL[icomp]=0.4;
    EXH[icomp]=0.4;
  }
  gPar = new TGraphAsymmErrors(nbinstot,X,Y,EXL,EXH,EYL,EYH);
  gPar->SetFillColor(6);
  gPar->SetLineColor(6);
  gPar->SetMarkerStyle(8);
  gPar->Draw("a2");
  gPar->Draw("p");
  return;
}



void histoCompare::showFitDiff(int isamp,int ibin,int iatt){

  // get (normalized) histogram with parameter modifications
  hMod = (TH1D*)hManager->getSumHistogramMod(isamp,ibin,iatt)->Clone("hmod");

  // draw MC histograms
  hMod->SetLineColor(kBlue);
  // subtract data from modified expectation
  hMod->Add(hManager->hData[isamp][ibin][iatt],-1.);
  hMod->Draw("h");
 
  //
  return;
}



void histoCompare::showFitResult(int isamp,int ibin,int iatt){

  // x axis lables
  TString xtitle[NATTMAX];
  xtitle[0] = "fiTQun e/#mu PID Par.";
  xtitle[1] = "fiTQun e/#pi^{0} PID Par.";
  xtitle[2] = "fiTQun mu/#pi PID Par.";
  xtitle[3] = "fiTQun Ring-Counting Par.";

  // xaxis title size
  double xtitlesize = 0.05;
  double ytitlesize = 0.05;
  double xlabelsize = 0.045;
  double ylabelsize = 0.045;

  // get (normalized) histogram with parameter modifications
  hMod = (TH1D*)hManager->getSumHistogramMod(isamp,ibin,iatt)->Clone("hmod");

  // get nominal histgram (also normalized)
  hTmp = hManager->getSumHistogram(isamp,ibin,iatt,1);

  // setup labels
  hMod->GetXaxis()->SetTitle(xtitle[iatt].Data());
  hMod->GetXaxis()->SetTitleSize(xtitlesize);
  hMod->GetXaxis()->SetLabelSize(xlabelsize);

  // draw MC histograms
  hMod->SetLineColor(kBlue);
  hTmp->SetLineColor(kRed);
 
  // draw data histograms
  double xmin = hManager->hData[isamp][ibin][iatt]->GetBinLowEdge(hManager->nBinBuffer+1);
  int nbinstot = hManager->hData[isamp][ibin][iatt]->GetNbinsX();
  double xmax =  hManager->hData[isamp][ibin][iatt]->GetBinLowEdge(nbinstot-hManager->nBinBuffer);
  xmax += hManager->hData[isamp][ibin][iatt]->GetBinWidth(nbinstot);
  hManager->hData[isamp][ibin][iatt]->GetXaxis()->SetRangeUser(xmin,xmax);
  hManager->hData[isamp][ibin][iatt]->SetMarkerStyle(8);
  hManager->hData[isamp][ibin][iatt]->SetTitle(Form("Detector Region %d Sample %d",ibin,isamp));
  hManager->hData[isamp][ibin][iatt]->GetXaxis()->SetTitle(xtitle[iatt].Data());
  hManager->hData[isamp][ibin][iatt]->GetXaxis()->SetTitleSize(xtitlesize);
  hManager->hData[isamp][ibin][iatt]->GetYaxis()->SetTitleSize(ytitlesize);
  hManager->hData[isamp][ibin][iatt]->GetXaxis()->SetLabelSize(xlabelsize);
  hManager->hData[isamp][ibin][iatt]->GetYaxis()->SetLabelSize(ylabelsize);
  hManager->hData[isamp][ibin][iatt]->GetXaxis()->CenterTitle(1);
  hManager->hData[isamp][ibin][iatt]->GetYaxis()->CenterTitle(1);

  hManager->hData[isamp][ibin][iatt]->Draw("e");

  hTmp->SetTitle(Form("DR_%d",ibin+1));
  hTmp->Draw("samehist");
  hMod->SetTitle(Form("DR_%d",ibin+1));
  hMod->Draw("samehist");

  //
  return;
}



//Show the effect of varying a single set of bias and smear parameters
void histoCompare::showFitEffect(int isamp,int ibin,int icomp,int iatt){

  // get the modified version of the specified histogram
  hMod = hManager->getModHistogram(isamp,ibin,icomp,iatt); //gets the modified histogram

  // get the default sum of histograms
  hTmp = hManager->getSumHistogram(isamp,ibin,iatt); //get the sum histogram

  // add default histograms to modified histograms
  for (int jcomp=0;jcomp<nComp;jcomp++){
    if (jcomp!=icomp){
      hMod->Add(hManager->hMC[isamp][ibin][jcomp][iatt]);
    }
  }

  // modified histogram is blue
  hMod->SetLineColor(kBlue);

  // default is red
  hTmp->SetLineColor(kRed);

  // scale to data
  double thenorm = hManager->normFactor;
  hMod->Scale(thenorm);
//  hTmp->Scale(thenorm);

  // draw to same canvas
  hMod->Draw("h");
  hTmp->Draw("sameh");
  hManager->hData[isamp][ibin][iatt]->SetMarkerStyle(8);
  hManager->hData[isamp][ibin][iatt]->Draw("samee");
  return;
}



///////////////////////////////////////////////////////////////////////////////////
//Plots the specified histogram before and after modifications
void histoCompare::showFitHisto(int isamp,int ibin,int icomp,int iatt){
  double smear = thePars->getHistoParameter(ibin,icomp,iatt,0);
  double bias  = thePars->getHistoParameter(ibin,icomp,iatt,1);
  cout<<"SMEAR: "<<smear<<endl;
  cout<<"BIAS:  "<<bias<<endl;
  hMod = hManager->getModHistogramMC(isamp,ibin,icomp,iatt);
  hMod->SetLineColor(kBlue);
  hMod->Draw("eh");
  hTmp = hManager->getHistogram(isamp,ibin,icomp,iatt);
  hTmp->SetLineColor(kRed);
  hTmp->Draw("sameeh");
  return;
}



void histoCompare::lnLWrapper(int& ndim, double* gout, double& result, double par[], int flg){

  for (int ipar=0;ipar<staticthis->thePars->nTotPars;ipar++){
    staticthis->thePars->setParameter(ipar,par[ipar]);
  }

  result = (double)staticthis->getTotLnL();
}


////////////////////////////////////////////
// Pre-fit parameters using MINUIT
void histoCompare::LnLPreFit(){

  //setup static this so wrapper doesn't segfault
  staticthis = this;

  //threshold to determine if a parameter is fit or not
  //if the MC histograms for this parameter have a size less than this value,
  //don't bother fitting them!
  double nthresh = 100.;

  //sets the precision of the fits
  double parerr = 0.1;  
  
  //individually fit each parameter
  int parindex =0;

  //parameter name container
  TString parnametmp;  

  //total number of parameters to be fit
  int npars = thePars->nTotPars;

  cout<<"$$$$$$$$$$$$$$$ LNL PRE FIT $$$$$$$$$$$$$$$"<<endl;
  cout<<"  ---------------------------------------- "<<endl;
  cout<<"  NUMBER OF PARAMETERS: "<<npars<<endl;
  cout<<"  PRECISION:            "<<parerr<<endl;
  cout<<"  ---------------------------------------  "<<endl;

  //fix parameters with too few events to be fit
  parindex = 0;
  for (int ibin=0;ibin<nBin;ibin++){
    for (int iatt=0;iatt<nAtt;iatt++){
      for (int icomp=0;icomp<nComp;icomp++){
        //name parameters
        parnametmp = Form("par_%d_",parindex);
        parnametmp.Append(binName[ibin].Data());
        parnametmp.Append("_");
        parnametmp.Append(compName[icomp].Data());
        parnametmp.Append("_");
        parnametmp.Append(attName[iatt].Data());
        parnametmp.Append("_");
        parnametmp.Append("smear");
        parName[ibin][icomp][iatt][0]=parnametmp.Data();
        parnametmp = Form("par_%d_",(parindex+1));
        parnametmp.Append(binName[ibin].Data());
        parnametmp.Append("_");
        parnametmp.Append(compName[icomp].Data());
        parnametmp.Append("_");
        parnametmp.Append(attName[iatt].Data());
        parnametmp.Append("_");
        parnametmp.Append("bias");
        parName[ibin][icomp][iatt][1]=parnametmp.Data();
        //get summed histogram
        hTot = (TH1D*)hManager->hMC[0][ibin][icomp][iatt]->Clone("htot");
        for (int isamp=1;isamp<nSamp;isamp++){
          hTot->Add(hManager->hMC[isamp][ibin][icomp][iatt]);
        }
        cout<<"total entries: "<<hTot->GetEntries()<<endl;
        //check to make sure histogram is above fitting threshold
        if (hTot->GetEntries()<nthresh){
          cout<<"  FIXING PARAMETER:  "<<parName[ibin][icomp][iatt][0].Data()<<" (# ENTRIES TOO LOW!) "<<endl; 
          fixPar[ibin][icomp][iatt][0]=1;
          thePars->fixParameter(ibin,icomp,iatt,0);
          thePars->fixParameter(ibin,icomp,iatt,1);
          cout<<"  FIXING PARAMETER:  "<<parName[ibin][icomp][iatt][1].Data()<<" (# ENTRIES TOO LOW!) "<<endl; 
        }
        hTot->Delete();
        parindex+=2;
      }
    }
  }
  cout<<"  ----------------------------------------  "<<endl;

  //setup the fitter!
  TFitter* fit = new TFitter(npars);

  //shut fitter up
  {
    double pp = -1;
    fit->ExecuteCommand("SET PRINTOUT",&pp,1);
  }

  //specify function to be fit
  fit->SetFCN(lnLWrapper);

  //set parameters to inital values
  TString pname;
  for (int ipar=0;ipar<thePars->nTotPars;ipar++){
    pname = Form("parameter%d",ipar);
    fit->SetParameter(ipar,pname.Data(),thePars->pars[ipar],parerr,0,0);
  }

  //fix all parameters
  for (int jpar=0;jpar<npars;jpar++){
    fit->FixParameter(jpar);
  }

  //release and fit normalization parameters
  for (int isyspar=(thePars->nTotPars-thePars->nNormPars);isyspar<thePars->nTotPars;isyspar++){
    if ((thePars->fixPar[isyspar])!=1)fit->ReleaseParameter(isyspar);
  }
  fit->ExecuteCommand("SIMPLEX",0,0);


  //fix all parameters
  for (int jpar=0;jpar<npars;jpar++){
    fit->FixParameter(jpar);
  }

  //release and fit flux and xsec parameters
  parindex = thePars->nTotPars-thePars->nSysPars;
  for (int isyspar=0;isyspar<(thePars->nSysPars-thePars->nNormPars);isyspar++){
    if ((thePars->fixPar[parindex])!=1)fit->ReleaseParameter(parindex);
    parindex++;
  }
  fit->ExecuteCommand("SIMPLEX",0,0);

  //fix all parameters
  for (int jpar=0;jpar<npars;jpar++){
    fit->FixParameter(jpar);
  }

  //run individual bias fits     
  parindex=0;
  for (int jbin=0;jbin<nBin;jbin++){
    for (int jatt=0;jatt<nAtt;jatt++){
      for (int jcomp=0;jcomp<nComp;jcomp++){
          //release biasparameter to be fit
          if (thePars->checkFixFlg(jbin,jcomp,jatt,1)==1){
            continue;
          }
          // get index of the smear parameter
          parindex = thePars->getParIndex(jbin,jcomp,jatt,1);
          fit->ReleaseParameter(parindex);
          cout<<"fitting "<<jbin<<jcomp<<jatt<<1<<" # "<<thePars->getParIndex(jbin,jcomp,jatt,1)<<endl;
          fit->ExecuteCommand("SIMPLEX",0,0);
          fit->FixParameter(parindex);
      }
    }
  }
  parindex=0;

  //run individual smear fits
  for (int jbin=0;jbin<nBin;jbin++){
    for (int jatt=0;jatt<nAtt;jatt++){
      for (int jcomp=0;jcomp<nComp;jcomp++){
          if (thePars->checkFixFlg(jbin,jcomp,jatt,0)==1){
            continue;
          }
          // get index of the smear parameter
          parindex = thePars->getParIndex(jbin,jcomp,jatt,0);
          fit->ReleaseParameter(parindex);
          cout<<"fitting "<<jbin<<jcomp<<jatt<<0<<" # "<<thePars->getParIndex(jbin,jcomp,jatt,0)<<endl;
          fit->ExecuteCommand("SIMPLEX",0,0);
          fit->FixParameter(parindex);
      }
    }
  }

  //fix all parameters
  for (int jpar=0;jpar<npars;jpar++){
    fit->FixParameter(jpar);
  }
  
  
  //print final results
  for (int ipar=0;ipar<npars;ipar++){
    thePars->setParameter(ipar,fit->GetParameter(ipar));
    cout<<"  PAR "<<ipar<<" FIT RESULT: "<<thePars->pars[ipar]<<endl;
  }
  cout<<"$$$$$$$$$$$$$ END LNL PRE FIT $$$$$$$$$$$$$"<<endl;

  return;
}


void histoCompare::sysParFit(){
  //setup static this so wrapper doesn't segfault
  staticthis = this;
  //sets the precision of the fits
  double parerr = 0.001;  
  cout<<"$$$$$$$$$$$$$$$$$ LNL SINLGE PAR FIT $$$$$$$$$$$$$$$$$"<<endl;
  cout<<"  ---------------------------------------- "<<endl;
  cout<<"  PRECISION:            "<<parerr<<endl;
  cout<<"  ---------------------------------------  "<<endl;

  //setup the fitter!

  //total number of parameters to be fit
  int npars = thePars->nTotPars;
  TFitter* fit = new TFitter(npars);
  //shut fitter up
  {
    double pp = 0;
    fit->ExecuteCommand("SET PRINTOUT",&pp,1);
  }
  //specify function to be fit
  fit->SetFCN(lnLWrapper);

  //setup parameters
  TString aname;
  for (int ipar=0;ipar<npars;ipar++){
    aname = "parameter";
    aname.Append(Form("_%d",ipar));
    fit->SetParameter(ipar,aname.Data(),thePars->pars[ipar],parerr,0,0);
  }
 
  //fix all params
  for (int jpar=0;jpar<npars;jpar++){
    fit->FixParameter(jpar);
  }
  //release single parameter to fit
  for (int jpar=(thePars->nTotPars-thePars->nSysPars);jpar<thePars->nTotPars;jpar++){
    fit->ReleaseParameter(jpar);
  }
  //fit that thang
  fit->ExecuteCommand("SIMPLEX",0,0); 
  fit->ExecuteCommand("SIMPLEX",0,0); 
  //print results
  for (int jpar=(thePars->nTotPars-thePars->nSysPars);jpar<thePars->nTotPars;jpar++){
    cout<<"PAR: "<<jpar<<" "<<fit->GetParameter(jpar)<<endl;
  }

  return;
}



/////////////////////////////////////////////////
// Fits only parameter "ipar"
void histoCompare::singleParFit(int ipar){
  //setup static this so wrapper doesn't segfault
  staticthis = this;
  //sets the precision of the fits
  double parerr = 0.001;  
  cout<<"$$$$$$$$$$$$$$$$$ LNL SINLGE PAR FIT $$$$$$$$$$$$$$$$$"<<endl;
  cout<<"  ---------------------------------------- "<<endl;
  cout<<"  PARAMETER: "<<ipar<<endl;
  cout<<"  PRECISION:            "<<parerr<<endl;
  cout<<"  ---------------------------------------  "<<endl;

  //setup the fitter!

  //total number of parameters to be fit
  int npars = thePars->nTotPars;
  TFitter* fit = new TFitter(npars);
  //shut fitter up
  {
    double pp = 0;
    fit->ExecuteCommand("SET PRINTOUT",&pp,1);
  }
  //specify function to be fit
  fit->SetFCN(lnLWrapper);

  //setup parameters
  TString aname;
  for (int kpar=0;kpar<npars;kpar++){
    aname = "parameter";
    aname.Append(Form("_%d",kpar));
    fit->SetParameter(kpar,aname.Data(),thePars->pars[kpar],parerr,0,0);
  }
 
  //fix all params
  for (int jpar=0;jpar<npars;jpar++){
    fit->FixParameter(jpar);
  }
  //release single parameter to fit
  fit->ReleaseParameter(ipar);
  //fit that thang
  fit->ExecuteCommand("SIMPLEX",0,0); 
  fit->ExecuteCommand("SIMPLEX",0,0); 
  //print results
  cout<<"RESULT: "<<fit->GetParameter(ipar)<<endl;
 
  return;
}


void histoCompare::printFitResults(const char* directory){

  cc = new TCanvas("cc","cc",700,900);
  cc->Divide(2,3);
  for (int isamp=0;isamp<nSamp;isamp++){
    for (int iatt=0;iatt<nAtt;iatt++){
      for (int ibin=0;ibin<nBin;ibin++){
        cc->cd(ibin+1);
        showFitResult(isamp,ibin,iatt);
        TString plotname = directory;
        plotname.Append(hTmp->GetName());
        plotname.Append(".png");
        cc->Print(plotname.Data()); 
      }
      TString plotname = directory;
      plotname.Append(Form("plt_samp%d_att%d",isamp,iatt));
      plotname.Append(".png");
      cc->Print(plotname.Data()); 
    }
  }

  return;
}



/////////////////////////////////////////////////////////////
// Fits all parameters using MINUIT
void histoCompare::LnLFit(){
  //setup static this so wrapper doesn't segfault
  staticthis = this;

  //sets the precision of the fits
  double parerr = 0.01;  
  
  //individually fit each parameter
  int parindex =0;

  //parameter name container
  TString parnametmp;  

  //will we fix all the smearing parameters?
  if (flgFixAllSmearPars){
    thePars->fixAllSmearPars();
  }

  //total number of parameters to be fit
  int npars = thePars->nTotPars;

  /////////////////////////////////////////////////////////
  //run the prefit
  LnLPreFit();

  ////////////////////////////////////////////////////////
  //setup the fitter!
  cout<<"$$$$$$$$$$$$$$$$$ LNL FIT $$$$$$$$$$$$$$$$$"<<endl;
  cout<<"  ---------------------------------------- "<<endl;
  cout<<"  NUMBER OF PARAMETERS: "<<npars<<endl;
  cout<<"  PRECISION:            "<<parerr<<endl;
  cout<<"  ---------------------------------------  "<<endl;

  TFitter* fit = new TFitter(npars);
  //shut fitter up
  {
    double pp = -1;
    fit->ExecuteCommand("SET PRINTOUT",&pp,1);
  }
  
  //specify function to be fit
  fit->SetFCN(lnLWrapper);
  //setup parameters
  TString aname;
  for (int ipar=0;ipar<npars;ipar++){
    aname = "parameter_";
    aname.Append(Form("%d",ipar));
    fit->SetParameter(ipar,aname.Data(),thePars->pars[ipar],parerr,0,0);
  } 
  parindex = 0;

  ///////////////////////////////////////////////////////////////////
  //do individual fits
  for (int jbin=0;jbin<nBin;jbin++){
    for (int jatt=0;jatt<nAtt;jatt++){

      //start of fit block//
      //
      //fix all parameters
      for (int jpar=0;jpar<npars;jpar++){
        fit->FixParameter(jpar);
      }

      //for bin 6, fit the xsec and flux parameters as well//
      if ((jbin==5)&&(jatt==0)){
        //release all flux and xsec pars
        for (int isyspar=(thePars->nTotPars-thePars->nSysPars);isyspar<thePars->nTotPars;isyspar++){
          if ((thePars->fixPar[isyspar])!=1)fit->ReleaseParameter(isyspar);  
          cout<<"fitting parameter: "<<isyspar<<endl;
         // fit->ReleaseParameter(isyspar);
        }
        //fit these pars first
        fit->ExecuteCommand("SIMPLEX",0,0);
        cout<<"--------------------------------"<<endl;
      }

      //release bias parameters
      for (int jcomp=0;jcomp<nComp;jcomp++){
        if (thePars->checkFixFlg(jbin,jcomp,jatt,1)!=1){
          parindex = thePars->getParIndex(jbin,jcomp,jatt,1);
          fit->ReleaseParameter(thePars->getParIndex(jbin,jcomp,jatt,1));
          cout<<"fitting parameter: "<<parindex<<endl;
        }
      }
      //run fit for bias parameters
      fit->ExecuteCommand("SIMPLEX",0,0);
      cout<<"--------------------------------"<<endl;

      //now release smear parameters
      for (int jcomp=0;jcomp<nComp;jcomp++){
        if (thePars->checkFixFlg(jbin,jcomp,jatt,0)!=1){
          parindex = thePars->getParIndex(jbin,jcomp,jatt,0);
          fit->ReleaseParameter(thePars->getParIndex(jbin,jcomp,jatt,0));
          cout<<"fitting parameter: "<<parindex<<endl;          
        }
      }

      fit->ExecuteCommand("SIMPLEX",0,0); //run the fit for ALL parameters
      cout<<"--------------------------------"<<endl;
      //end of fit block//
      
    }
  }

  //set final results
  for (int ipar=0;ipar<npars;ipar++){
    thePars->setParameter(ipar,fit->GetParameter(ipar));
  }

  //calculate rough errors and print results
  calcRoughParErr();

  //end
  cout<<"  ----------------------------------------  "<<endl;
  cout<<"$$$$$$$$        FIT  COMPLLETE      $$$$$$$$"<<endl;

  return;


}



void histoCompare::getTotLnL1D(double& result,int npar, double par[]){

  for (int ipar=0;ipar<npar;ipar++){
    thePars->setParameter(ipar,par[ipar]);
  }

  result = getTotLnL();
}



///////////////////////////////////////////////////
// Compute the likelihood component coming from priors
double histoCompare::getPriorLnL(){

  double priorLnL = 0.;
  double pull = 0;
  for (int isys=0;isys<thePars->nSysPars;isys++){
    pull = thePars->getSysParameter(isys)-1.;
    pull/=thePars->sysParUnc[isys];
    priorLnL+=(0.5)*(pull*pull);
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  //Contribution from bias and smear priors
  // If using gaussion priors on bias and smear parameters, evelauate the likelihood here.
  // This is done by calling an atmfit pars method that will sum the contributions 
  if (flgUsePriorsInFit){
    double parpriorlnl = thePars->calcLogPriors();
    priorLnL += parpriorlnl; 
#ifdef VERBOSE
    cout<<"prior LnL: "<<parpriorlnl<<endl;
    cout<<"Total LnL: "<<totL<<endl;
#endif

  }


  return priorLnL;
}



////////////////////////////////////////////////
//Compute the total log liklihood by comparing all histograms
///*
double histoCompare::getTotLnL(){

  double totL = 0.;

  nDOF = 0.;
  ////////////////////////////////////////
  //contribution from histogram comparison  
  for (int isamp=0;isamp<nSamp;isamp++){
//  for (int isamp=0;isamp<1;isamp++){
    for (int ibin=0;ibin<nBin;ibin++){
      for (int iatt=0;iatt<nAtt;iatt++){
       	TH1D* hPrediction = (TH1D*)hManager->getSumHistogramMod(isamp,ibin,iatt,1); //< get normalized histogram.
      	TH1D* hDataTmp = (TH1D*)hManager->getHistogramData(isamp,ibin,iatt);
        double partialLnL = hManager->histoLogL;
#ifdef VERBOSE
        cout<<"partial LnL for sample: "<<isamp<<" bin: "<< ibin << " att: "<<iatt<<endl;
        cout<<"partial LnL: "<<partialLnL<<endl;
        cout<<"total LnL: "<<totL<<endl;
#endif
	      totL+=partialLnL;
        nDOF+=hManager->nDOF;
      }
    }
  }

  //////////////////////////////////////////////
  //contribution from flux/xsec priors
  double pull;
#ifndef T2K
  double priorlnl = getPriorLnL();
#ifdef VERBOSE
  cout<<"prior LnL: "<<priorlnl<<endl;
#endif
  totL+=priorlnl;
  //cout<<"prior: "<<priorlnl<<endl;
//  for (int isys=0;isys<thePars->nSysPars;isys++){
//    pull = thePars->getSysParameter(isys)-1.;
//    pull/=thePars->sysParUnc[isys];
//    totL+=(0.5)*(pull*pull);
//  }

#else
  for (int isys=0;isys < 2;isys++){
    pull = thePars->sysPar[isys]-1.;
    pull/=thePars->sysParUnc[isys];
    totL+=(0.5)*(pull*pull);
  }
  totL += thePars->cov->getLikelihood();
#endif




  return totL;
  
}
//*/


double histoCompare::getTotSumSq(){
  double totsumsq = 0.;
  return totsumsq;
}



double histoCompare::getSumSq(TH1D* h1, TH1D* h2){
  double sumsq = 0.;
  double diff;
  for (int ibin=10;ibin<=(h1->GetNbinsX()-10);ibin++){
    diff = h1->GetBinContent(ibin)-h2->GetBinContent(ibin);
    sumsq += (diff*diff);
  }
  return sumsq;
}



////////////////////////////////////////////////////////////////////
//evalute log-likelihood between two histograms
double histoCompare::getLnL(TH1D* h1, TH1D* h2){
  double lnL = 0.;
//  double diff;
//  double term;
  double c1; //data
  double c2; //mc
  double errmc; //mcerr
  //double norm = hManager->normFactor; //normalization
  double norm = 1.0;
//  double dof=0.;
//  double quaderr;



  ///////////////////////////////////////////////////
  //assume poisson errors
 // for (int ibin=10;ibin<=(10);ibin++){
  int nedgebins= 5; //< need to ignore some edge bins
  for (int ibin=nedgebins;ibin<=(h1->GetNbinsX()-nedgebins);ibin++){
    c1 = h1->GetBinContent(ibin); //MC
    c2 = h2->GetBinContent(ibin); //data
    errmc = h1->GetBinError(ibin);

    if (c2<7) continue;
     
    lnL+=evalLnL(c2,c1,norm); //< tn186 likelihood definition
  //  lnL+=evalGausChi2WithError(c2,c1,errmc); //< tn186 likelihood definition

}
  return lnL;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
// Read in histograms from a specified file
void histoCompare::readFromFile(const char* filerootname,int nsamp, int nbin, int ncomp, int natt){

  // set variables
  nSamp = nsamp;
  nBin  = nbin;
  nComp = ncomp; 
  nAtt  = natt;

  // read in histograms
  hManager = new histoManager(filerootname,nsamp,nbin,ncomp,natt);
  
  //
  return;
}


void histoCompare::readFromFile(const char* filerootname,int nsamp, int nbin, int ncomp, int nmode, int natt){
  nSamp = nsamp;
  nBin  = nbin;
  nComp = ncomp; 
  nAtt  = natt;
  nMode = nmode;
  hManager = new histoManager(filerootname,nsamp,nbin,ncomp,natt,nmode,true);
  double ndataevents=0;
  double nmcevents=0;
  double events;
}


void histoCompare::makeFakeData(){
   // sets all modified histograms as data 
   for (int isamp=0; isamp<nSamp; isamp++){
     for (int ibin=0; ibin<nBin; ibin++){
       for (int iatt=0; iatt<nAtt; iatt++){
         TString hname = hManager->hData[isamp][ibin][iatt]->GetName();
         hManager->hData[isamp][ibin][iatt] = (TH1D*)hManager->getSumHistogramMod(isamp,ibin,iatt,0)
         ->Clone(hname.Data());
       }
     }
   }
   hManager->normFactor = 1.0;
   return;
}

/////////////////////////////////////////////////////////////////////
//initializes all necessary compoents
void  histoCompare::initialize(histoManager* hm, atmFitPars* apars){
  cout<<"histoCompare: Initialization: "<<endl;
  thePars = apars;
  hManager = hm;
  nComp = thePars->nComponents;
  nBin  = thePars->nBins;
  nAtt  = thePars->nAttributes;
  nSamp = thePars->nSamples;
  cout<<"    bins: "<<nBin<<endl;
  cout<<"    sampless: "<<nSamp<<endl;
  cout<<"    components: "<<nComp<<endl;
  cout<<"    attributes: "<<nAtt<<endl; 
  return; 
}

void histoCompare::setupPars(int nsyspars){
  thePars = new atmFitPars(nSamp,nBin,nComp,nAtt,"tn186");
 // thePars->setNorm(Norm);
  hManager->setFitPars(thePars);
  return;
}

void histoCompare::setupPars(atmFitPars *a)
{
  thePars = a;
  hManager->setFitPars(thePars);
}

histoCompare::histoCompare(){
  nameTag = "manualSetup";
  nMCHist=0;
  nDataHist=0;
  cout<<"created comparison object: "<<nameTag.Data()<<endl;
  //setup initial debug
  int jhist = 0;
  cc = new TCanvas("cc","cc",700,600);
  while (jhist<10){
      parDebug[jhist][0] = 1.0;
      parDebug[jhist][1] = 0.0;
      jhist++;
  }
  useLnLType=0;
  return;
}


//////////////////////////////////////////////////////////
// Use this constructor 
histoCompare::histoCompare(const char* parfile, bool sep)
  : separateNeutMode(sep)
{

  //nominal MCMC output file
  MCMCOutputFile = "mcmctree.root";

  //read in parameter file
  runPars = new sharedPars(parfile);
  runPars->readParsFromFile();
  nameTag= runPars->globalRootName.Data();

  //setup canvas
  cc = new TCanvas("cc","cc",700,600);
  
  //set LnL definition to default (0)
  useLnLType=0;

  //MCMC tuning parameter
  tunePar = runPars->MCMCTunePar;

  //Use smear parameters or no?
  flgFixAllSmearPars = runPars->flgFixAllSmearPars;

  //MCMC nsteps;
  MCMCNSteps = runPars->MCMCNSteps;

  //MCMC nburn
  MCMCNBurnIn = runPars->MCMCNBurnIn;

  //read in pre-filled histograms using histoManager
  int nbins = runPars->nFVBins;
  int ncomponents = runPars->nComponents;
  int nsamples = runPars->nSamples;
  int nattributes  = runPars->nAttributes;
  TString histofilename = runPars->hFactoryOutput;
  readFromFile(histofilename.Data(),nsamples,nbins,ncomponents,nattributes); 
  hManager->nBinBuffer = runPars->nBinBuffer;

  //setup fit parameters
  thePars = new atmFitPars(parfile);
  hManager->setFitPars(thePars);
  int nsyspars = thePars->nSysPars;
  if (flgFixAllSmearPars){
    thePars->fixAllSmearPars();
  }

  
  //read in splines if you're into that
  if (runPars->useSplinesFlg){
    setupSplines(runPars->splineFactoryOutput.Data());
  };

  //Should we use priors in this fit?
  flgUsePriorsInFit = runPars->flgUsePriorsInFit;
  // if we're using priors, get them
  if (flgUsePriorsInFit){
    for (int iatt=0; iatt<nattributes; iatt++){
      double smearwidth =  runPars->kr->getKeyD(Form("smearPriorWidthAtt%d",iatt));
      double biaswidth =  runPars->kr->getKeyD(Form("biasPriorWidthAtt%d",iatt));
      cout<<"setting prior for attribute "<<iatt<<"to: "<<smearwidth<<" "<<biaswidth<<endl;
      thePars->setHistoParPrior(iatt,0,smearwidth);
      thePars->setHistoParPrior(iatt,1,biaswidth);
    }
  }

  //Read in parameters from previous fit?
  if (runPars->flgUseFitParFile){
    thePars->readPars(runPars->kr->getKeyS("fitParFile"));
  }

  //Are there physical bounds?
  if (runPars->flgUsePhysLoBound){
    cout<<"histoCompare: Using physical bounds!"<<endl;
    for (int iatt=0; iatt<nattributes; iatt++){
      double physbound = runPars->kr->getKeyD(Form("physLoBoundAtt%d",iatt));
      if (physbound>=0){
        hManager->setLoBound(iatt,physbound);      
      }
    }
    return;
  }

}




void histoCompare::makeResidualErrorMaps(const char* outdir){

  // get directory name
  TString dirname = outdir;

  // make a new file
  TString filename = dirname.Data();
  filename.Append("ResidualErrs.root");
  TFile *resfile = new TFile(filename.Data(),"RECREATE");

  // for histo names
  TString hname;
  TH2FV* hfv[10];
  int hindex=0;

  // make histograms (this may have to be adjusted depending
  // on the cut variables used)
  
 
  // NuE PID
  hfv[hindex] = new TH2FV("ResUncNuEPID",0);
  int attindex = 0;
  int sampleindex = 0;
  double cutval = 0.;
  int accept_greater = 1;
  for (int ibin=1; ibin<=hfv[0]->GetNumberOfBins(); ibin++){
  
   // get normalized MC
   TH1D* hmc = hManager->getSumHistogramMod(sampleindex,ibin-1,attindex,1);

   // get data
   TH1D* hdata = hManager->hData[sampleindex][ibin-1][attindex];

   // calc residual
   double res = calcResErrCumFracDiff(hmc,hdata,cutval,accept_greater);
//   double res = calcFracDiff(hmc,hdata);

   // set bin content
   hfv[hindex]->SetBinContent(ibin,res);

  }
  hfv[hindex]->Write();
  hindex++;

  // NuE PI0
  hfv[hindex] = new TH2FV("ResUncNuEPi0",0);
  attindex = 1;
  sampleindex = 0;
  cutval = 0.;
  accept_greater = 0;
  for (int ibin=1; ibin<=hfv[0]->GetNumberOfBins(); ibin++){
  
   // get normalized MC
   TH1D* hmc = hManager->getSumHistogramMod(sampleindex,ibin-1,attindex,1);

   // get data
   TH1D* hdata = hManager->hData[sampleindex][ibin-1][attindex];

   // calc residual
   double res = calcResErrCumFracDiff(hmc,hdata,cutval,accept_greater);

   // set bin content
   hfv[hindex]->SetBinContent(ibin,res);

  }
  hfv[hindex]->Write();
  hindex++;

  // NuMu PID
  hfv[hindex] = new TH2FV("ResUncNuMuPID",0);
  attindex = 0;
  sampleindex = 1;
  cutval = 0.;
  accept_greater = 0;
  for (int ibin=1; ibin<=hfv[hindex]->GetNumberOfBins(); ibin++){
  
   // get normalized MC
   TH1D* hmc = hManager->getSumHistogramMod(sampleindex,ibin-1,attindex,1);

   // get data
   TH1D* hdata = hManager->hData[sampleindex][ibin-1][attindex];

   // calc residual
   double res = calcResErrCumFracDiff(hmc,hdata,cutval,accept_greater);

   // set bin content
   hfv[hindex]->SetBinContent(ibin,res);

  }
  hfv[hindex]->Write();
  hindex++;

  //
  resfile->Close();

  //
  return;

}









#endif
