#include "histoCompare.h"
#include <stdlib.h>
using namespace std;

////////////////////////////////////////////
// Run a markov chian
//   Usage: runMCMCFit <parfile1> <randomseed2> <outputfilename3>
int main(int argc, char* argv[]){
 
// cout<<"helloworld"<<endl;

 // setup random
 int rand_seed = atoi(argv[2]);
 randy = new TRandom2(rand_seed);
 cout<<"TRandom seed is: "<<randy->GetSeed()<<endl;

 // get par file
 TString parfilename = argv[1];
 cout<<"Parameter file "<<parfilename.Data()<<endl;

 // get ouput file (4th argument)
 TString outfilename = argv[3];

 // make new histogram comparison object
 histoCompare* hc= new histoCompare(parfilename.Data());
 
 // set up the fit
 hc->MCMCOutputFile = outfilename.Data();
// hc->thePars->fixAll
// hc->hManager->setLoBound(3,0);
// hc->hManager->setLoBound(2,0);
// hc->thePars->readPars("/nfs/data41/t2k/amissert/atmos/head/atmFitTools/pars/lnlfit_priors.root");
// hc->thePars->fixAllAttPars(3);
 // run the mcmc
 hc->runMCMC(-1); 

 //
 return 0;
 

}



