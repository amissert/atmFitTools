#include "eventSelectors.h"
#include "calcEnu.h"
#include "mcLargeArray.h"
#include "TH1D.h"
#include "TMath.h"
#include "TChain.h"
#include "TString.h"


using namespace std;


class testLargeArray{
  public:

  testLargeArray(const char* dir);
  TChain* chmc;
  mcLargeArray* fastevents;
  TH1D* henu_numu;
  TH1D* henu_nue;

  void fillhistos();
};


testLargeArray::testLargeArray(const char* dir){
  
  chmc = new TChain("h1");
  chmc->Add(dir);
  fastevents = new mcLargeArray(chmc,10000000);
}

void testLargeArray::fillhistos(){
  int nbins = 100;
  double xmin = 0;
  double xmax = 2000;
  henu_numu = new TH1D("henu_numu","henu_numu",nbins,xmin,xmax);
  henu_nue = new TH1D("henu_nue","henu_nue",nbins,xmin,xmax);
  for (int iev=0; iev<chmc->GetEntries(); iev++){
    int ipass = selectNuMu(fastevents->vnhitac[iev],
                fastevents->vfqnsubev[iev],
                fastevents->vfqenumu[iev],
                fastevents->vfqemom[iev],
                fastevents->vfqmumom[iev],
                fastevents->vfqpid[iev],
                fastevents->vfqnring[iev]);

    if (!ipass) continue;    
    if (fastevents->vpassnumu[iev]){
      henu_numu->Fill(fastevents->vfqenumu[iev],fastevents->vweight[iev]);
      henu_nue->Fill(fastevents->vfqenue[iev],fastevents->vweight[iev]);
    }
  }

  return;
}

