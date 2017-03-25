#ifndef UNIFORMLIST_H
#define UNIFORMLIST_H

#include "TRandom2.h"
#include "TMath.h"
#include <vector>
#include <algorithm>

#define LISTMAX 1000000


using namespace std;



class uniformList{

 public:

 // constructor
 uniformList(int npts, int nmin, int nmax){
   nPts = npts;
   nMin = nmin;
   nMax = nmax;
   init();
 }

 int  nMax;
 int  nMin;
 int  nPts;
 int  nSeed;
 vector<int> vlist;
 int  getAt(int ipoint){return vlist.at(ipoint);}
 void init(){
   float dx =  ((float)nMax - (float)nMin)
               /((float)nPts-1);
   float x=(float)nMin;
   for (int i=0; i<nPts; i++){
     int point = (int)(x+0.5);
     vlist.push_back(point);
     x+=dx;
   } 
   return;
 }

};


#endif

