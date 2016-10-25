#include "histoCompare.h"
#include <stdlib.h>
using namespace std;


int main(int argc, char* argv[]){

 // setup random
 int rand_seed = atoi(argv[2]);
// randy->SetSeed(rand_seed);
 randy = new TRandom2(rand_seed);
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
 hc->hManager->setLoBound(3,0);
 hc->hManager->setLoBound(2,0);
 hc->thePars->readPars("/nfs/data41/t2k/amissert/atmos/head/atmFitTools/pars/initial_lnlfit_pars.root");
 hc->setPar(26,1.0);
 hc->setPar(27,0.0);
 hc->setPar(154,1.0);
 hc->setPar(155,0.0);
 hc->setPar(106,1.0);
 hc->setPar(107,0.0);
 hc->thePars->fixPar[26];
 hc->thePars->fixPar[27];
 hc->thePars->fixPar[154];
 hc->thePars->fixPar[155];
 hc->thePars->fixPar[106];
 hc->thePars->fixPar[107];
 hc->tunePar = 0.055;


 // run the mcmc
 hc->runMCMC(-1); 

 //
 return 0;

}



