#ifndef PREPROCESS_C
#define PREPROCESS_C

#include "preProcess.h"

/////////////////////////////////////////////////////////////////
// Setup the FV bin histogram for getFVBin()
void preProcess::setFVBinHisto(){
  hFVBins = new TH2FV("hfvbins",1);
}


/////////////////////////////////////////////////////////////////
//Create a TGraph that is used to make weights to see effect of
//correcting for cosmic flux mis-modeling
void  preProcess::setWeightHistogram(const char* file, const char * name){

  //read in histogram
  TFile* hfile = new TFile(file);
  hWeight = (TH1D*)hfile->Get(name);
 
  //convert to graph
  const int nn = hWeight->GetNbinsX();
  double xx[nn];
  double yy[nn];
  for (int i=0;i<nn;i++){
    xx[i] = hWeight->GetBinCenter(i+1);
    yy[i] = hWeight->GetBinContent(i+1);
  }
  gWeight = new TGraph(nn,xx,yy);

  gWeight->Draw("alc");

  useWeights = 1;

  return;
}

////////////////////////////////////////////////////////////
//Takes in a chain and loops over all files in the chain
//For each file, a new file with a modified TTree is created
void preProcess::processAllFiles(TChain* chain){
  int nfiles = chain->GetNtrees();
  TObjArray* listOfFiles = chain->GetListOfFiles();
  TString tag;
  TString fname;
  for (int ifile=0;ifile<nfiles;ifile++){
    tag = Form("_%d_",ifile);
    fname = listOfFiles->At(ifile)->GetTitle();
    processFile(fname,tag);
  }
  return;
}

////////////////////////////////////////
//sets up pointers given a TChain
void preProcess::setTree(TChain* chin){
  TTree* trin = (TTree*)chin;
  setTree(trin);
  return;
}

///////////////////////////////////////
//sets up pointers
void preProcess::setTree(TTree* trin){
  tr = trin;
  fq = new fqEvent(tr);
  vis = new visRing(fq); 
  return;
}

///////////////////////////////////////////////////
// reads in a file and processes the h1 tree inside
void preProcess::processFile(const char* fname,const char* outname){

  //make output file name
  TString outputName = outDir.Data(); //name of directory
  outputName.Append(nameTag.Data()); //global name
  outputName.Append(outname); //number of this file
  outputName.Append(".root");

  //get existing h1 tree
  cout<<"  opening file: "<<fname<<endl;
  TFile* fin = new TFile(fname);
  TTree* intree = (TTree*)fin->Get("h1");
  cout<<"  got tree: "<<intree->GetEntries()<<endl;
  setTree(intree); //< set pointers to current tree

  //make new tree
  cout<<"  create file: "<<outputName.Data()<<endl;
  fout = new TFile(outputName.Data(),"recreate");
  setupNewTree(); 

  //fill new tree
  int neventsnew = preProcessIt();
  cout<<"  filled it with "<<neventsnew<<" events!" <<endl;


  //clean up
  if (trout->GetEntries()>0) fout->Write();
  fin->Close();
  fout->Close();
  
  return;   
}

///////////////////////////////////
//Gets a weight for an event
//Usefull for making fake data sets
float preProcess::getWeight(){
//  absmode = TMath::Abs(fq->mode);
//  float enu   = fq->pmomv[0];
  evtweight = 1.0;
  //CCQE norm bin1 
//  if ((absmode==1)&&(enu<200.)){
//    evtweight = 1.5;
//  }
  //CCQE norm bin2 
//  if ((absmode==1)&&(enu>200.)&&(enu<400.)) evtweight*=1.2;
  //CCQE norm bin3 
//  if ((absmode==1)&&(enu>400.)&&(enu<800.)) evtweight*=0.9;
  //CCQE norm bin4 
//  if ((absmode==1)&&(enu>800.)) evtweight*=1.05;


  if (useWeights){
    evtweight = gWeight->Eval(fq->fq1rmom[0][2],0,"s");
  }

  return evtweight;
}



void preProcess::setParFileName(const char* fname){
  
  // set file name
  parFileName=fname;
  
  //read in parameters!
  cout<<"preProcess: Reading parameters from file "<<parFileName.Data()<<endl;
  sharedPars* runpars = new sharedPars(parFileName.Data());
  runpars->readParsFromFile();
  nameTag = runpars->globalRootName;
  cout<<"nametag: "<<nameTag.Data()<<endl;
  FVBinning = runpars->preProcessFVBinning; //< flag for FV binning type in getBin()
  setFVBinHisto();
  MCComponents = runpars->preProcessMCComponents; //< flag for MC component definitions in getComponent()
  MCSamples = runpars->preProcessMCSamples; //< flag for MC sample definitions in getSample()
  NHITACMax = runpars->preProcFCCut;
  EVisMin = runpars->preProcEVisCut;
  WallMin = runpars->preProcWallMinCut;
  ToWallMin = runpars->preProcToWallMinCut;
  NSEMax = runpars->preProcNseMax0;
  NSEMin = runpars->preProcNseMin;
  InGateMin = runpars->preProcInGateCut; 
  flgAddMoreVars = runpars->preProcAddMoreVars;

  // list of attributes to use
  nAttributes = runpars->nAttributes;
  attributeList[0] = runpars->fQAttName0;
  attributeList[1] = runpars->fQAttName1;
  attributeList[2] = runpars->fQAttName2;
  attributeList[3] = runpars->fQAttName3;
  attributeList[4] = runpars->fQAttName4;
  attributeList[5] = runpars->fQAttName5;
  attributeList[6] = runpars->fQAttName6;
  attributeList[7] = runpars->fQAttName7;


}



///////////////////////////////////////////////
//calculates the FV bin for an event
int preProcess::getBin(){

  ////////////////////////////////////////////////////////////
  //calculate fiducial volume variables
  //use electron hypothesis
  TVector3 vpos;
  vpos.SetXYZ(fq->fq1rpos[0][2][0],fq->fq1rpos[0][2][1],fq->fq1rpos[0][2][2]);
  TVector3 vdir;
  vdir.SetXYZ(fq->fq1rdir[0][2][0],fq->fq1rdir[0][2][1],fq->fq1rdir[0][2][2]);
  wall = calcWall2(&vpos);
  towall = calcToWall(&vpos,&vdir);
  // calculate additional fv variables as well
  for (int isubev=0; isubev<fq->fqnse; isubev++){
    vpos.SetXYZ(fq->fq1rpos[isubev][2][0],fq->fq1rpos[isubev][2][1],fq->fq1rpos[isubev][2][2]);
    vdir.SetXYZ(fq->fq1rdir[isubev][2][0],fq->fq1rdir[isubev][2][1],fq->fq1rdir[isubev][2][2]);
    fq1rwall[isubev][2] = calcWall2(&vpos);
    fq1rtowall[isubev][2] = calcToWall(&vpos,&vdir);
    vpos.SetXYZ(fq->fq1rpos[isubev][1][0],fq->fq1rpos[isubev][1][1],fq->fq1rpos[isubev][1][2]);
    vdir.SetXYZ(fq->fq1rdir[isubev][1][0],fq->fq1rdir[isubev][1][1],fq->fq1rdir[isubev][1][2]);
    fq1rwall[isubev][1] = calcWall2(&vpos);
    fq1rtowall[isubev][1] = calcToWall(&vpos,&vdir);
  }
  //true towall
  for (int ipart=0; ipart<fq->npar; ipart++){
    vpos.SetXYZ(fq->posv[0],fq->posv[1],fq->posv[2]);
    vdir.SetXYZ(fq->dirv[ipart][0],fq->dirv[ipart][1],fq->dirv[ipart][2]);
    towallv[ipart]=calcToWall(&vpos,&vdir);
    wallv2=calcWall2(&vpos);
  }
  // even more FV related variables
  if (flgAddMoreVars>0){
     // add 1r perimiter and mincone	  
     for (int isubev=0; isubev<fq->fqnse; isubev++){
       vpos.SetXYZ(fq->fq1rpos[isubev][2][0],fq->fq1rpos[isubev][2][1],fq->fq1rpos[isubev][2][2]);
       vdir.SetXYZ(fq->fq1rdir[isubev][2][0],fq->fq1rdir[isubev][2][1],fq->fq1rdir[isubev][2][2]);
       fq1rperim[isubev][2] = calcPerimeter(&vpos,&vdir);
       fq1rmincone[isubev][2] = calcMinCone(&vpos,&vdir);
       vpos.SetXYZ(fq->fq1rpos[isubev][1][0],fq->fq1rpos[isubev][1][1],fq->fq1rpos[isubev][1][2]);
       vdir.SetXYZ(fq->fq1rdir[isubev][1][0],fq->fq1rdir[isubev][1][1],fq->fq1rdir[isubev][1][2]);
       fq1rperim[isubev][1] = calcPerimeter(&vpos,&vdir);
       fq1rmincone[isubev][1] = calcMinCone(&vpos,&vdir);
    }
    //true perimeter and mincone
    if (flgAddMoreVars>1){
      for (int ipart=0; ipart<fq->npar; ipart++){
        vpos.SetXYZ(fq->posv[0],fq->posv[1],fq->posv[2]);
        vdir.SetXYZ(fq->dirv[ipart][0],fq->dirv[ipart][1],fq->dirv[ipart][2]);
        perimv[ipart]=calcPerimeter(&vpos,&vdir);
        minconev[ipart]=calcMinCone(&vpos,&vdir);
      }
    }
  }

  ////////////////////////////////
  // separate into bins by FVBinning parameter

  //////////////////////////////////////////
  //"simple" FV Binning for atm
  if (FVBinning==0){
    if ((wall<200.)&&(wall>80.)) return 1;
    if (wall<80.) return 0;
    return 2;
  }

  ////////////////////////////////////////
  // cosmic binning by entering surface
  if (FVBinning==1){
    double Rrec = TMath::Sqrt(pow(fq->fq1rpos[0][2][0],2)+pow(fq->fq1rpos[0][2][1],2));
    double Zrec = fq->fq1rpos[0][2][2];
    double Zcut = 1410;
    double Rcut = 1290;
    if ((Zrec>Zcut)&&(Rrec<Rcut)) return 0; //< top entering
    if ((Zrec<Zcut)) return 1; //< side entering
    if ((Zrec>Zcut)&&(Rrec>Rcut)) return 2; //< corner entering
  }

  //////////////////////////////////
  // all in one bin
  if (FVBinning==2){return 0;}
  
  //////////////////////////////////
  // towall binning for cosmics
  if (FVBinning==3){
    if (towall<500.) return 0;
    if ((towall>-500)&&(towall<1000)) return 1;
    if (towall>=1000) return 2;
  }

  //////////////////////////////////////
  // binning using wall/towallhistogram
  if (FVBinning==4){
      int fvbin = hFVBins->FindBin(towall,wall)-1;
      return fvbin;
   }
//    if (fvbin<0){
//      cout<<"Bad FV value:"<<endl;
//      cout<<"  wall:   "<<wall<<endl;
//      cout<<"  towall: "<<towall<<endl;
//  //  }
      //return fvbin;
 //   if (fq->fq1rmom[0][1]<230.) return fvbin;
  //  if (fq->fq1rmom[0][1]>800.) return fvbin+12;
   // return fvbin+6;
  

  return -1;
}


/////////////////////////////////
//Simple initial cuts
int preProcess::passCuts(){

  /////////////////////
  //tmp cuts
 // if (towallv[0]<80.) return 0;

  /////////////////////
  //Fully Contained Cut
  if ((int)fq->nhitac>NHITACMax) return 0;

  //////////////////////
  //Visible Energy Cut
  if (fq->fq1rmom[0][1]<EVisMin) return 0;

  ////////////////
  //FV Basic Cuts
  if (wall<WallMin) return 0; 
  if (towall<ToWallMin) return 0;  

  /////////////////////////
  //Number of subevent cuts
  if (fq->fqnse>NSEMax) return 0;
  if (fq->fqnse<NSEMin) return 0;
 
  /////////////
  //In-gate cut
  if (InGateMin>0){
    double tdecay = fq->fq1rt0[1][1]-fq->fq1rt0[0][2];
    if (tdecay<InGateMin) return 0;
  }

  ////////////////////////
  // optional masking cut
  if (flgUseSpikeMask>0){
     if (!passMask(hmask,fq1rwall[0][2])) return 0;
  }

  ////////////////////
  //all cuts passed!!
  return 1;
}

///////////////////////////////
//returns the # of subevents-1
int preProcess::getSample(){

  //atmospheric selections
  if (MCSamples==0){
    if (fq->fqnse==1) return 0;
    if (fq->fqnse==2) return 1;
    if (fq->fqnse>2)  return 2;
  }
  
  
  //cosmic selection
  if (MCComponents==1){
//    double Rrec = TMath::Sqrt(pow(fq->fq1rpos[0][2][0],2)+pow(fq->fq1rpos[0][2][1],2));
//    double Zrec = fq->fq1rpos[0][2][2];
//    double Zcut = 1410;
//    double Rcut = 1290;
//    if ((Zrec>Zcut)&&(Rrec<Rcut)) return 0; //< top entering
//    if ((Zrec<Zcut)) return 1; //< side entering
//    if ((Zrec>Zcut)&&(Rrec>Rcut)) return 2; //< corner entering
    return 0;
  }

  //hybrid pi0 selection
  if (MCSamples==2){
//    double Rrec = TMath::Sqrt(pow(fq->fq1rpos[0][2][0],2)+pow(fq->fq1rpos[0][2][1],2));
//    double Zrec = fq->fq1rpos[0][2][2];
//    double Zcut = 1410;
//    double Rcut = 1290;
//    if ((Zrec>Zcut)&&(Rrec<Rcut)) return 0; //< top entering
//    if ((Zrec<Zcut)) return 1; //< side entering
//    if ((Zrec>Zcut)&&(Rrec>Rcut)) return 2; //< corner entering
    return 0;
  }

  cout<<"preProcess:  Error sample not defined!"<<endl;
  return -1;
}

//////////////////////////////////////
//get code for MC true component type
int preProcess::getComponent(){

  ////////////////////////////
  // useful for cuts
  absmode = TMath::Abs(fq->mode);
  int absnu   = TMath::Abs(fq->ipnu[0]);
 
  /////////////////////////////////////////
  // visible + NEUT event selection for atm
  if (MCComponents==0){
    if ((absmode>0)&&(absmode<30)){
      if ((vis->nve==1)&&((vis->nvp==0)&&(vis->nvmu==0)&&(vis->nvpi0==0)&&(vis->nvpip==0))) return 0; //CC single e
      if ((vis->nvmu==1)&&((vis->nvp==0)&&(vis->nvpi0==0)&&(vis->nvpip==0))) return 1; //CC single mu
      if ((vis->nve==1)&&(vis->nvmu==0)) return 2; //CC e + other
      if ((vis->nve==0)&&(vis->nvmu==1)) return 3; //CC mu + other
      return 4; //CC Other
    }
    else{
      if ((vis->nvpi0>0)) return 5;  //single pi0
      return 6;  //NC other
    }
  }

  /////////////////////////////////////////
  // visible only components for atm
  if (MCComponents==2){

    // single electron
    if ((vis->nve==1)&&(vis->nvp==0)&&(vis->nvmu==0)&&(vis->nvpi0==0)&&(vis->nvpip==0)) return 0;
    // single muon
    if ((vis->nve==0)&&(vis->nvp==0)&&(vis->nvmu==1)&&(vis->nvpi0==0)&&(vis->nvpip==0)) return 1;
    // electron + other
    if (vis->nve==1) return 2;
    // muon + other
    if (vis->nvmu==1) return 3;
    // pi0 with no other
    if ((vis->nvpi0==1)&&(vis->nvpip==0)) return 4;
    // other
    return 5;
  }
 
  //////////////////////////////////////////
  // cosmic selectoin
  if (MCComponents==1){
    return 0;
  }
  

  //////////////////////////////////////////
  // hybrid pi0 selectoin
  if (MCComponents==3){
    return 0;
  }

  cout<<"preProcess: ERROR MC component not defined!";
  return -1;
}

//loop over all events and sort into bins, samples and components
int preProcess::preProcessIt(){
  int nev = tr->GetEntries();
  int naccepted = 0;
  for (int i=0;i<nev;i++){
    //get info for event
    if ((i%1000)==0) cout<<"event:  "<<i<<endl;
    tr->GetEntry(i);
    //calc FV bin and fill FV variables
    nbin=getBin();
    if (nbin<0.) continue;
    //apply cuts
    if (!passCuts()) continue;
    naccepted++;
    vis->fillVisVar(); //get visible ring information
    fillAttributes(fq);
    ncomponent=getComponent();
    nsample=getSample();
    evtweight=getWeight();
    trout->Fill();
  }

  return naccepted;
}

///////////////////////////////////////
//returns the index of the best 2R fit
int preProcess::getBest2RFitID(){
  
  int nfits = fq->fqnmrfit;

  double ngLnLBest = 10000000.;
  int bestindex = 0;

  for (int ifit=0;ifit<nfits;ifit++){
    int fitID = TMath::Abs(fq->fqmrifit[ifit]); //< fit fit ID code
//    int diff = (TMath::Abs(fitID-20000000));
//    cout<<"diff: "<<diff<<endl;
//    cout<<"ifit: "<<ifit<<endl;
    if ((TMath::Abs(fitID-20000000))>100) continue; //< we want best 2R fits
    if (fq->fqmrnll[ifit]<ngLnLBest){
      bestindex = ifit;
      ngLnLBest=fq->fqmrnll[ifit];
    }
  }
  best2RID = fq->fqmrifit[bestindex];
  return bestindex;
}

////////////////////////////////////////
//fills fiTQun attribute array
void preProcess::fillAttributes(fqEvent* fqevent){

  // Fill the cmap that matches attribute names to values
  fillAttributeMap(fqevent);

  // fill the attribute[] array with the values you want to use for this analysis
  for (int i=0; i<nAttributes; i++){
    double value = attributeMap[attributeList[i].Data()];
    attribute[i] = value;
  }


/*  attribute[0] = fqevent->fq1rnll[0][2]-fqevent->fq1rnll[0][1];
  int ibest = getBest2RFitID();
  double best1Rnglnl = fmin(fqevent->fq1rnll[0][1],fqevent->fq1rnll[0][2]);
  attribute[1] = best1Rnglnl-fqevent->fqmrnll[ibest];
  attribute[2] = wall;
  attribute[5] = fqevent->fq1rnll[1][2]-fqevent->fq1rnll[1][1];
//  attribute[2] = fqevent->fq1rmom[0][2];
  attribute[3] = fqevent->fq1rpos[0][2][2];
  double xx = fqevent->fq1rpos[0][2][0];
  double yy = fqevent->fq1rpos[0][2][1];
  attribute[4] = TMath::Sqrt(xx*xx + yy*yy);
  */

  return;
}


////////////////////////////////////////
//fills cmap of possible fitqun attributes
void preProcess::fillAttributeMap(fqEvent* fqevent){

  // PID e vs. mu ratio of first subevent
  attributeMap["fqelike"] = fqevent->fq1rnll[0][2]-fqevent->fq1rnll[0][1];

  // Ring Counting (RC) parameter
  int ibest = getBest2RFitID();
  double best1Rnglnl = fmin(fqevent->fq1rnll[0][1],fqevent->fq1rnll[0][2]);
  fqrcpar = best1Rnglnl-fqevent->fqmrnll[ibest];
  attributeMap["fqrcpar"] = fqrcpar;

  // Reconstructed distance from wall
  attributeMap["fqwall"] = wall;

  // Reconstructed momentum (muon)
  attributeMap["fq1rmumom"] = fqevent->fq1rmom[0][2];

  // Reconstructed momentum (electron)
  attributeMap["fq1remom"] = fqevent->fq1rmom[0][1];

  // pi0 likelihood
  attributeMap["fqpi0like"] = fqevent->fq1rnll[0][1] - fqevent->fqpi0nll[0];

  // pi0 mass
  attributeMap["fqpi0mass"] = fqevent->fqpi0mass[0];

  // pi0 photon angle
  attributeMap["fqpi0photangle"] = fqevent->fqpi0photangle[0];

  // pi0 wall min
  TVector3 vpos;
  vpos.SetXYZ(fqevent->fqpi0pos[0][0],fqevent->fqpi0pos[0][1],fqevent->fqpi0pos[0][2]);
  TVector3 vdir1;
  vdir1.SetXYZ(fqevent->fqpi0dir1[0][0],fqevent->fqpi0dir1[0][1],fqevent->fqpi0dir1[0][2]);
  TVector3 vdir2;
  vdir2.SetXYZ(fqevent->fqpi0dir2[0][0],fqevent->fqpi0dir2[0][1],fqevent->fqpi0dir2[0][2]);
  double pi0wall= calcWall2(&vpos);
  double pi0towall1 = calcToWall(&vpos,&vdir1);
  double pi0towall2 = calcToWall(&vpos,&vdir2);
  attributeMap["fqpi0wall"] = pi0wall;

  // pi0 towall min
  attributeMap["fqpi0towallmin"] = fmin(pi0towall1,pi0towall2);

  // pi0 total momentum
  attributeMap["fqpi0momtot"] = fqevent->fqpi0momtot[0];

  return;
}

void preProcess::setupNewTree(){
  tr->SetBranchStatus("*",0);
  tr->SetBranchStatus("fq*",1);
  tr->SetBranchStatus("*v",1);
  tr->SetBranchStatus("ipnu",1);
  //tr->SetBranchStatus("cluster*",1);
  tr->SetBranchStatus("mode",1);
  tr->SetBranchStatus("nhitac",1);
  trout = tr->CloneTree(0); //clone but don't copy data
  trout->CopyAddresses(tr); //set addresses
  
  //add new branches
  trout->Branch("attribute",attribute,"attribute[1000]/F");
  trout->Branch("fqrcpar",&fqrcpar,"fqrcpar/F");
  trout->Branch("ncomponent",&ncomponent,"ncomponent/I");
  trout->Branch("nsample",&nsample,"nsample/I");
  trout->Branch("nbin",&nbin,"nbin/I");
  trout->Branch("nvis",&vis->nvis,"nvis/I");
  trout->Branch("nvmu",&vis->nvmu,"nvmu/I");
  trout->Branch("nve",&vis->nve,"nve/I");
  trout->Branch("nvgam",&vis->nvgam,"nvgam/I");
  trout->Branch("nvpip",&vis->nvpip,"nvpip/I");
  trout->Branch("nvpi0",&vis->nvpi0,"nvpi0/I");
  trout->Branch("nvp",&vis->nvp,"nvp/I");
  trout->Branch("nvk",&vis->nvk,"nvk/I");
  trout->Branch("nvoth",&vis->nvoth,"nvoth/I");
  trout->Branch("vispid",vis->vispid,"vispid[100]/I");
  trout->Branch("fqwall",&wall,"fqwall/F");
  trout->Branch("fqtowall",&towall,"fqtowall/F");
  trout->Branch("fq1rwall",fq1rwall,"fq1rwall[10][7]/F");
  trout->Branch("fq1rtowall",fq1rtowall,"fq1rtowall[10][7]/F");
  trout->Branch("towallv",towallv,"towallv[50]");
  trout->Branch("wallv2",&wallv2,"wallv2");
  trout->Branch("evtweight",&evtweight,"evtweight/F");
  trout->Branch("best2RID",&best2RID,"best2RID/I");
  trout->Branch("fq1rperim",fq1rperim,"fq1rperim[10][7]/F");
  trout->Branch("fq1rmincone",fq1rmincone,"fq1rmincone[10][7]/F");

  return;
}


/////////////////////////////
//empty constructor
preProcess::preProcess(){
  nFiles=0;
  useWeights=0;
}


//////////////////////////////////////////
//read in parameters and run preprocessing!
void preProcess::runPreProcessing(){
  
  //read in parameters!
  cout<<"preProcess: Reading parameters from file "<<parFileName.Data()<<endl;
  sharedPars* runpars = new sharedPars(parFileName.Data());
  runpars->readParsFromFile();
  nameTag = runpars->globalRootName;
  cout<<"nametag: "<<nameTag.Data()<<endl;
  FVBinning = runpars->preProcessFVBinning; //< flag for FV binning type in getBin()
  if (FVBinning==4) setFVBinHisto();
  MCComponents = runpars->preProcessMCComponents; //< flag for MC component definitions in getComponent()
  MCSamples = runpars->preProcessMCSamples; //< flag for MC sample definitions in getSample()
  NHITACMax = runpars->preProcFCCut;
  EVisMin = runpars->preProcEVisCut;
  WallMin = runpars->preProcWallMinCut;
  ToWallMin = runpars->preProcToWallMinCut;
  NSEMax = runpars->preProcNseMax0;
  NSEMin = runpars->preProcNseMin;
  InGateMin = runpars->preProcInGateCut; 
  flgAddMoreVars = runpars->preProcAddMoreVars;
  flgUseSpikeMask = runpars->preProcMaskFlg;
  if (flgUseSpikeMask>0){
     TString fname = runpars->preProcMaskFile.Data();
     TFile* maskfile = new TFile(fname.Data());
     cout<<"preProc: Getting spike mask from file: "<<fname.Data()<<endl;     
     hmask = (TH1D*)maskfile->Get("hmask");
  }

  // list of attributes to use
  nAttributes = runpars->nAttributes;
  attributeList[0] = runpars->fQAttName0;
  attributeList[1] = runpars->fQAttName1;
  attributeList[2] = runpars->fQAttName2;
  attributeList[3] = runpars->fQAttName3;
  attributeList[4] = runpars->fQAttName4;
  attributeList[5] = runpars->fQAttName5;
  attributeList[6] = runpars->fQAttName6;
  attributeList[7] = runpars->fQAttName7;


  //create data and mc chains
  chmc = new TChain("h1");
  chdat = new TChain("h1");
  cout<<"preProc: adding MC files: "<<runpars->preProcessFilesMC.Data()<<endl;
  cout<<"preProc: adding Data files: "<<runpars->preProcessFilesData.Data()<<endl;
  chmc->Add(runpars->preProcessFilesMC.Data());
  chdat->Add(runpars->preProcessFilesData.Data());
  if (chmc->GetEntries()<1){
    cout<<"preProcess WARNING: no events in MC chain"<<endl;
  }
  if (chdat->GetEntries()<1){
    cout<<"preProcess WARNING: no events in Data chain"<<endl;
  }

  //process the files
  outDir = runpars->preProcessOutDir.Data();
  nameTag.Append("_ppmc");
  processAllFiles(chmc);
  nameTag = runpars->globalRootName.Data();
  nameTag.Append("_ppdata");
  processAllFiles(chdat); 

  cout<<"preProcess: Complete!"<<endl;

  //////////////////////////
  return;
}

#endif
