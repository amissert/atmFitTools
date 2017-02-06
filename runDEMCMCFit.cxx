#include "histoCompare.h"
#include <stdlib.h>

/////////////////////////////////////////////////////////////////////////////////////
// Runs a differential evoloution markov chian
//
// Usage:
//    runDEMCMCFIT <parameterfile1> <differentialfiles2> <randomseed3> <outputfile4>
using namespace std;


int main(int argc, char* argv[]){

 // setup random
 int rand_seed = atoi(argv[3]);
 randy->SetSeed(rand_seed);
 cout<<"TRandom seed is: "<<randy->GetSeed()<<endl;

 // get par file
 TString parfilename = argv[1];
 cout<<"Parameter file "<<parfilename.Data()<<endl;

 // get ouput file (4th argument)
 TString outfilename = argv[4];
 cout<<"Output file "<<outfilename.Data()<<endl;

 // make new histogram comparison object
 histoCompare* hc= new histoCompare(parfilename.Data());

 // names of files of mcmc point differences
 TString diffnames = argv[2]; 
 diffnames.ReplaceAll("\\","");
 cout<<"Differential files "<<diffnames.Data()<<endl;

 // set up the fit
 hc->MCMCOutputFile = outfilename.Data();
 cout<<"hc output file: "<<hc->MCMCOutputFile.Data()<<endl;
// hc->hManager->setLoBound(2,0);
// hc->hManager->setLoBound(3,0);
// hc->thePars->readPars("/nfs/data41/t2k/amissert/atmos/head/atmFitTools/pars/initial_lnlfit_pars.root");
 hc->diffChainFileName = diffnames.Data();
// hc->tunePar = 0.02;

 // fix these pars
// hc->thePars->fixAllAttPars(3);

 // run the mcmc
 hc->runDEMCMC(-1); 

 //
 return 0;

}



