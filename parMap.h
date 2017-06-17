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

  // get name
  TString getSystParName(int ipar);

  private:

  void init();
};



#endif  
