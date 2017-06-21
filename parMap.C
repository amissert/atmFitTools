#ifndef  parMap_C
#define  parMap_C

#include "parMap.h"


int parMap::getSysIndex(int isys){
  int istart = nParsTot - nSyst;
  return istart + isys;
}

int parMap::getGlobalIndex(int ibin, int icomp, int iatt, int itype){
  return gIndex[ibin][icomp][iatt][itype];
}


TString parMap::getName(int ipar){

  if (ipar > nParsTot){
    cout<<"Error: only "<<nSyst<<" total parameters"<<endl;
  }

  return pName[ipar];
}


TString parMap::getSystParName(int ipar){

  if (ipar > nSyst){
    cout<<"Error: only "<<nSyst<<" systematic parameters"<<endl;
  }

  return systName[ipar];
}


parMap::parMap(int nbin, int ncomp, int natt, int nsyst){

  nSyst = nsyst;
  nBin = nbin;
  nComp = ncomp;
  nAtt = natt;
  nParsTot = nbin*ncomp*natt*2 + nsyst;
  cout<<"Total parameters:"<<nParsTot<<endl;

  init();
}


void parMap::init(){
  
  // done  by hand
  systName[0]="CCQE par. 1";    
  systName[1]="CCQE par. 2";    
  systName[2]="CCQE par. 3";    
  systName[3]="CCQE par. 4";    
  systName[4]="CCQE par. 5";    
  systName[5]="CCQE par. 6";    
  systName[6]="CCQE par. 7";    
  systName[7]="CCQE par. 8";    
  systName[8]="CCQE par. 9";    
  systName[9]="CCQE par. 10";    
  systName[10]="CCQE par. 11";    
  systName[11]="CCQE par. 12";    
  systName[12]="CCQE par. 13";    
  systName[13]="CCQE par. 14";    
  systName[14]="Sub-GeV flux";    
  systName[15]="Mulit-GeV flux";    
  systName[16]="CCnQE";    
  systName[17]="NC";    
  systName[18]="#mu/e ratio";    
  systName[19]="Region 0 Sample 0 norm.";
  systName[20]="Region 0 Sample 1 norm.";
  systName[21]="Region 0 Sample 2 norm.";
  systName[22]="Region 1 Sample 0 norm.";
  systName[23]="Region 1 Sample 1 norm.";
  systName[24]="Region 1 Sample 2 norm.";
  systName[25]="Region 2 Sample 0 norm.";
  systName[26]="Region 2 Sample 1 norm.";
  systName[27]="Region 2 Sample 2 norm.";
  systName[28]="Region 3 Sample 0 norm.";
  systName[29]="Region 3 Sample 1 norm.";
  systName[30]="Region 3 Sample 2 norm.";
  systName[31]="Region 4 Sample 0 norm.";
  systName[32]="Region 4 Sample 1 norm.";
  systName[33]="Region 4 Sample 2 norm.";
  systName[34]="Region 5 Sample 0 norm.";
  systName[35]="Region 5 Sample 1 norm.";
  systName[36]="Region 5 Sample 2 norm.";

  //
  TString compName[6];
  compName[0]="Single e ";
  compName[1]="Single #mu ";
  compName[2]="e + Other ";
  compName[3]="#mu + Other ";
  compName[4]="Single #pi^{0} ";
  compName[5]="Single Hadron ";

  //
  TString attName[6];
  attName[0] = "e/#mu PID ";
  attName[1] = "e/#pi^{0} PID ";
  attName[2] = "#mu/#pi PID ";
  attName[3] = "RC par. ";

  TString typeName[2];
  typeName[0]="#beta_{1}";
  typeName[1]="#beta_{0}";

  int parindex = 0;
  for (int ibin=0; ibin<nBin; ibin++){;
    for (int icomp=0; icomp<nComp; icomp++){
      for (int iatt=0; iatt<nAtt; iatt++){
        for (int itype=0; itype<2; itype++){
          gIndex[ibin][icomp][iatt][itype] = parindex;
//          TString name = Form("Region %d ",ibin);
          TString name = "";
          name.Append(compName[icomp].Data());
          name.Append(attName[iatt].Data());
          name.Append(typeName[itype].Data());
          cout<<"par #"<<parindex<<" gets "<<name.Data()<<endl;
          pName[parindex]=name.Data();
          parindex++;
        }
      }
    }    
  }

  for (int isys=0; isys<nSyst; isys++){
    pName[parindex] = systName[isys].Data();
    parindex++;
  }

  return;
}



#endif  

