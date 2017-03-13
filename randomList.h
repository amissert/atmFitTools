#ifndef RANDOMLIST_H
#define RANDOMLIST_H

#include "TRandom2.h"
#include "TMath.h"
#include <vector>
#include <algorithm>

#define LISTMAX 100000


using namespace std;



class randomList{

 public:

 // constructor
 randomList(int npts, int nmax, int nseed);

 int  nMax;
 int  nEvts;
 int  nSeed;
// int  eventnum[LISTMAX];
 vector<int> vlist;
 int  getAt(int ipoint);
 void init();

};

#ifdef CINTMODE
#include "randomList.cxx"
#endif

#endif

