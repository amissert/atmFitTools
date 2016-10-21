#include "histoCompare.h"
#include <stdlib.h>
using namespace std;


int main(int argc, char* argv[]){

 // setup random
 int rand_seed = atoi(argv[2]);
 randy->SetSeed(rand_seed);
// randy = new TRandom2(rand_seed);
 cout<<"TRandom seed is: "<<randy->GetSeed()<<endl;

 // get par file
 TString parfilename = argv[1];
 cout<<"Parameter file"<<parfilename.Data()<<endl;

 // get ouput file (4th argument)
 TString outfilename = argv[3];

 // make new histogram comparison object
 histoCompare* hc= new histoCompare(parfilename.Data());

 // set up the fit
 hc->MCMCOutputFile = outfilename.Data();
 hc->thePars->fixAllAttPars(1);
 hc->hManager->setLoBound(3,0);
 hc->hManager->setLoBound(4,0);
 hc->thePars->readPars("/nfs/data41/t2k/amissert/atmos/head/atmFitTools/initial_lnlfit.root");
 hc->tunePar = 0.044;


 // run the mcmc
 hc->runMCMC(-1); 

 //
 return 0;

}


/*
int rootrun(const char* parfilename){

 // make new histogram comparison object
 histoCompare* hc= new histoCompare(parfilename);

 // run the mcmc
 hc->runMCMC(-1); 

 //
 return 0;
}

*/



