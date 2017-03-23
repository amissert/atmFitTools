
#include "TString.h"
#include "TChain.h"
#include "TTree.h"
#include "TGraph.h"
#include "TMath.h"
#include "TRandom2.h"
#include "TVector3.h"
#include "TCanvas.h"
#include <iostream>
#include "../TH2FV.h"
#include "../FVCalculators.h"

const int RNDSEED = 1000;
TRandom2* randy = new TRandom2(RNDSEED);



///////////////////////////////////////////////////////////////////
double calcWallSphere(TVector3* vpos){

  double Rmax = 1690.;

  double rr = vpos->Mag();

  return (Rmax-rr);

}


///////////////////////////////////////////////////////////////////
double calcToWallSphere(TVector3* vpos, TVector3* vdir){

  double Rmax = 1690.;
  double r0 = vpos->Mag();
  if (r0>Rmax) return -10.;
  double Vr = vdir->Dot(*vpos);

  double root2 = - Vr + TMath::Sqrt( (Vr*Vr) - (r0*r0) + (Rmax*Rmax));
  if (root2>=0.) return root2;

  double root1 = - Vr - TMath::Sqrt( (Vr*Vr) - (r0*r0) + (Rmax*Rmax));
  if (root1>=0.) return root1; 


  cout<<"oops..."<<endl;

  return 0.;

}


///////////////////////////////////////////////////////////////////
TVector3* getRandomDir(){

  // sk geom
  double Rmax = 1690.;
  double Zmax = 1810.;

  // get random point on sphere
  double dir_theta = 2.*TMath::Pi()*randy->Rndm();
  double dir_phi = TMath::ACos(2.*randy->Rndm() - 1.) + (TMath::Pi()/2.);
 
  //get direction 
  TVector3* dir = new TVector3();
  dir->SetXYZ( TMath::Cos(dir_theta)*TMath::Cos(dir_phi),
              TMath::Sin(dir_theta)*TMath::Cos(dir_phi),
              TMath::Sin(dir_phi) );

  return dir;

}



///////////////////////////////////////////////////////////////////
TVector3* getRandomPos(){

  // sk geom
  double Rmax = 1690.;
  double Zmax = 1810.;
 
  // get random position
  double pos_z = (2.*Rmax*randy->Rndm() - Rmax); 
  double pos_r = Rmax*TMath::Sqrt(randy->Rndm());
  double pos_theta = 2.*TMath::Pi()*randy->Rndm();
  TVector3* pos = new TVector3();
  pos->SetXYZ( pos_r*TMath::Cos(pos_theta),
                pos_r*TMath::Sin(pos_theta),
                pos_z );

  return pos;

}


///////////////////////////////////////////////////////////////////
TH2D* fillPosHistoZ(int npoints){

  // sk geom
  double Rmax = 1690.;
  double Zmax = 1810.;

  TH2D* hrad = new TH2D("hrad","hrad",100,-Rmax,Rmax,100,-Zmax,Zmax);

  for (int ipt=0; ipt<npoints; ipt++){

    // get random position
    TVector3* pos = getRandomPos();
    
    // fill histogram
    hrad->Fill(pos->X(), pos->Z());

  }

  hrad->Draw("colz");

  return hrad;
}


///////////////////////////////////////////////////////////////////
TH2D* fillPosHistoR(int npoints){

  // sk geom
  double Rmax = 1690.;
  double Zmax = 1810.;

  TH2D* hrad = new TH2D("hrad","hrad",100,-Rmax,Rmax,100,-Rmax,Rmax);

  for (int ipt=0; ipt<npoints; ipt++){
    
    // get random position
    TVector3* pos = getRandomPos();
    
    // fill histogram
    hrad->Fill(pos->X(), pos->Y());

  }

  hrad->Draw("colz");

  return hrad;
}

///////////////////////////////////////////////////////////////////
TH2D* fillDirHistoXY(int npoints){

  TH2D* hrad = new TH2D("hrad","hrad",100,-1,1,100,-1,1);

  for (int ipt=0; ipt<npoints; ipt++){
    
    // get random position
    TVector3* dir = getRandomDir();
    
    // fill histogram
    hrad->Fill(dir->X(), dir->Y());

  }

  hrad->Draw("colz");

  return hrad;
}



///////////////////////////////////////////////////////////////////
TH2D* fillDirHistoXZ(int npoints){

  TH2D* hrad = new TH2D("hrad","hrad",100,-1,1,100,-1,1);

  for (int ipt=0; ipt<npoints; ipt++){
    
    // get random position
    TVector3* dir = getRandomDir();
    
    // fill histogram
    hrad->Fill(dir->X(), dir->Z());

  }

  hrad->Draw("colz");

  return hrad;
}



///////////////////////////////////////////////////////////////////
void fillRandomPoints(TH2FV* hfv, int npoints){

  for (int ipt=0; ipt<npoints; ipt++){
    
    // get random direction
    TVector3* dir = getRandomDir();

    // get random position
    TVector3* pos = getRandomPos();
   
    // get FV variables
    double wall = calcWall(pos);
    double towall = calcToWall(pos,dir);
   
    // fill histogram
    hfv->Fill(towall, wall);
    
    dir->Delete();
    pos->Delete();

  }

  hfv->Draw("colz");

  return;
}


///////////////////////////////////////////////////////////////////
void fillRandomPointsSphere(TH2FV* hfv, int npoints){

  for (int ipt=0; ipt<npoints; ipt++){
    
    // get random direction
    TVector3* dir = getRandomDir();

    // get random position
    TVector3* pos = getRandomPos();
   
    // get FV variables
    double wall = calcWallSphere(pos);
    double towall = calcToWallSphere(pos,dir);
   
    // fill histogram
    hfv->Fill(towall, wall);
    
    dir->Delete();
    pos->Delete();

  }

  hfv->Draw("colz");

  return;
}






