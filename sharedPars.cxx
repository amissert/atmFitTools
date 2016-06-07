<<<<<<< HEAD:sharedPars.cxx
#ifndef SHAREDPARS_C
#define SHAREDPARS_C

#include "TString.h"


#include "sharedPars.h"

double sharedPars::getParD(const char* parname){
  return kr->getKeyD(parname);
}

TString sharedPars::getParS(const char* parname){
  return kr->getKeyS(parname);
}

int sharedPars::getParI(const char* parname){
  return kr->getKeyI(parname);
}

void sharedPars::readParsFromFile(const char* filename){

  
  //create object to read in keys from file
  kr = new keyread(parFileName.Data());

  //read in contents
  kr->readFile();

  //set parameters to values
  useSplinesFlg = kr->getKeyI("useSplinesFlg");
  MCMCNSteps= kr->getKeyI("MCMCNSteps");
  MCMCTunePar = kr->getKeyD("MCMCTunePar");
  fixAllSmearFlg = kr->getKeyI("fixAllSmearFlg");
  FVBinName0= kr->getKeyS("FVBinName0");
  FVBinName1=kr->getKeyS("FVBinName1");
  FVBinName2=kr->getKeyS("FVBinName2");
  FVBinName3=kr->getKeyS("FVBinName3");
  FVBinName4=kr->getKeyS("FVBinName4");
  FVBinName5=kr->getKeyS("FVBinName5");
  fQAttName0=kr->getKeyS("fQAttName0");
  fQAttName1=kr->getKeyS("fQAttName1");
<<<<<<< HEAD:sharedPars.cxx
  fQAttName2=kr->getKeyS("fQAttName2");
  fQAttName3=kr->getKeyS("fQAttName3");
  fQAttName4=kr->getKeyS("fQAttName4");
  fQAttName5=kr->getKeyS("fQAttName5");
  fQAttName6=kr->getKeyS("fQAttName6");
  fQAttName7=kr->getKeyS("fQAttName7");
  fQAttName8=kr->getKeyS("fQAttName8");
=======
>>>>>>> xy_xy:sharedPars.cxx
  hFactoryDataFiles=kr->getKeyS("hFactoryDataFiles");
  hFactoryMCFiles=kr->getKeyS("hFactoryMCFiles");
  hFactoryOutput=kr->getKeyS("hFactoryOutput");
  MCComponentName0=kr->getKeyS("MCComponentName0");
  MCComponentName1=kr->getKeyS("MCComponentName1");
  MCComponentName2=kr->getKeyS("MCComponentName2");
  MCComponentName3=kr->getKeyS("MCComponentName3");
  MCComponentName4=kr->getKeyS("MCComponentName4");
  MCComponentName5=kr->getKeyS("MCComponentName5");
  MCComponentName6=kr->getKeyS("MCComponentName6");
  sampleName0=kr->getKeyS("sampleName0");
  sampleName1=kr->getKeyS("sampleName1");
  sampleName2=kr->getKeyS("sampleName2");
  sysParName0=kr->getKeyS("sysParName0");
  sysParName1=kr->getKeyS("sysParName1");
  sysParName2=kr->getKeyS("sysParName2");
  sysParName3=kr->getKeyS("sysParName3");
  sysParName4=kr->getKeyS("sysParName4");
  sysParName5=kr->getKeyS("sysParName5");
  sysParName6=kr->getKeyS("sysParName6");
  sysParName7=kr->getKeyS("sysParName7");
  sysParName8=kr->getKeyS("sysParName8");
  nFVBins     = kr->getKeyI("nFVBins");
  nSamples    = kr->getKeyI("nSamples");
  nComponents = kr->getKeyI("nComponents");
  nAttributes  = kr->getKeyI("nAttributes"); 
  nSysPars    = kr->getKeyI("nSysPars");
  preProcessFilesMC = kr->getKeyS("preProcessFilesMC"); 
  preProcessOutDir = kr->getKeyS("preProcessOutDir"); 
  preProcessFilesData = kr->getKeyS("preProcessFilesData"); 
  preProcessFilesBANFF = kr->getKeyS("preProcessFilesBANFF");
  preProcessFilesSpline = kr->getKeyS("preProcessFilesSpline");
  preProcessMCComponents = kr->getKeyI("preProcessMCComponents");
  preProcessFVBinning = kr->getKeyI("preProcessFVBinning");
  preProcessMCSamples = kr->getKeyI("preProcessMCSamples");
  preProcAddMoreVars = kr->getKeyI("preProcAddMoreVars");
  preProcMaskFile = kr->getKeyS("preProcMaskFile");
  preProcMaskFlg = kr->getKeyI("preProcMaskFlg");
  globalRootName = kr->getKeyS("globalRootName");
  splineFactoryOutput = kr->getKeyS("splineFactoryOutput");
  sysParType = kr->getKeyS("sysParType");
  NMCMCPts = kr->getKeyI("NMCMCPts");
  MCMCBurnIn=kr->getKeyI("MCMCBurnIn");
  NMCEvents=kr->getKeyI("NMCEvents");
  MCMCFile=kr->getKeyS("MCMCFile");
  preProcFCCut=kr->getKeyI("preProcFCCut");;
  preProcEVisCut=kr->getKeyD("preProcEVisCut");
  preProcWallMinCut=kr->getKeyD("preProcWallMinCut");
  preProcToWallMinCut=kr->getKeyD("preProcToWallMinCut");
  preProcNseMax0=kr->getKeyI("preProcNseMax");
  preProcNseMin=kr->getKeyI("preProcNseMin");
  preProcInGateCut=kr->getKeyD("preProcInGateCut");
  NDataEvents = kr->getKeyI("NDataEvents");
  flgUseNormPars = kr->getKeyI("flgUseNormPars");
}

sharedPars::sharedPars(const char* parfilename){
  parFileName = parfilename;
//  readParsFromFile();
}
