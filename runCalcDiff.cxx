#include "markovDiff.h"

#include <iostream>
#include <stdlib.h>

using namespace std;


////////////////////////////////////////////////////////////////////////////
// Loops over pairs of points in a MCMC cloud For each pair, find the vector
// difference between parameters and save to ouput file
//
// Usage:
//  ./runCalcMCMCDiff <mcmcfiles1> <nburn2> <nfill3> <randomseed4> <outputfile5>
//
int main(int argc, char* argv[]){

 // number of filled pairs
 int nfilledpairs = atoi(argv[3]);

 // setup random
 int rand_seed = atoi(argv[4]);
 randy->SetSeed(rand_seed);
 cout<<"TRandom seed is: "<<randy->GetSeed()<<endl;

 // get mcmc
 TString mcmcfiles = argv[1];
 cout<<"MCMC files are "<<mcmcfiles.Data()<<endl;

 // get ouput file (4th argument)
 TString outfilename = argv[5];

 // get nburn
 int nburnsteps = atoi(argv[2]);

 // make object
 markovDiff* md = new markovDiff(mcmcfiles,nburnsteps,1);
 md->outFileName = outfilename.Data(); //< set name of output file
 md->fillDiffPars(nfilledpairs);

 //
 return 0;

}
