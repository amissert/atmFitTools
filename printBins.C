
#include "TH2FV.cxx"
#include <iostream>

using namespace std;

void printBins(TH2FV* h){

  for (int i=0; i<h->GetNumberOfBins(); i++){
    double binc = h->GetBinContent(i);
    cout<<"bin #"<<i<<" has content "<<binc<<endl;
  }
}

