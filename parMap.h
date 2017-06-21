#ifndef  parMap_H
#define  parMap_H

#define NFITPAR_MAX 500

#include <iostream>
#include <map>
#include "TString.h"

using namespace std;



class parMap{

  public:

  // constructor
  parMap(int nbin=6, int ncomp=6, int natt=4, int nsyst=37);

  // data
  TString pName[NFITPAR_MAX]; 
  TString systName[NFITPAR_MAX]; 
  int nParsTot;
  int nSyst;
  int nBin;
  int nComp;
  int nAtt;
  int gIndex[6][6][4][2]; //< [nbins][ncomp][natt][ntype]

  // get name
  TString getSystParName(int ipar);
  TString getName(int ipar);
  int getGlobalIndex(int ibin, int icomp, int iatt, int itype);
  int getSysIndex(int isys);
  private:

  void init();
};



#endif  
