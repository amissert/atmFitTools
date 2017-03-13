#ifndef RANDOMLIST_CXX 
#define RANDOMLIST_CXX

#include "randomList.h"

using namespace std;

// constructor
randomList::randomList(int npts, int nmax, int nseed){
  
  nEvts = npts;
  nMax = nmax;
  nSeed = nseed;

  init();

  return;
}

// initialize
void randomList::init(){

  TRandom2* randgen = new TRandom2(nSeed);

  for (int ipt=0; ipt<=nEvts; ipt++){
    vlist.push_back(randgen->Integer(nMax));
  }

  // sort
  sort(vlist.begin(),vlist.end());

  return;
}

// git
int randomList::getAt(int ievt){
  return vlist.at(ievt);
}



#endif
