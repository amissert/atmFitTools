#ifndef CALCENU_CXX
#define CALCENU_CXX

#include "calcEnu.h"

float calcEnu(float lpmom, TVector3 ldir, int muflg){

 // t2k beam coordinates in SK frame
 float beam_x = 0.669764;
 float beam_y = -0.7421791;
 float beam_z = 0.024228;
 TVector3 beamdir;
 beamdir.SetXYZ(beam_x,beam_y,beam_z);

 // binding energy
 float Eb = 27.; //< in MeV 

 // masses
 float Me = 0.511;
 float Mmu = 105.658;
 float Mp = 938.28;
 float Mn = 939.57;

 // energy
 float Ml = 0.;
 if (muflg==1){
  Ml = Mmu;
 }
 else{
  Ml = Me;
 }

 float El = TMath::Sqrt( (lpmom*lpmom) + (Ml*Ml));
 // cosine to beam
// float angle = ldir.Angle(beamdir);
 float cosbeam = ldir.Unit().Dot(beamdir.Unit());
// float cosbeam = TMath::Cos(angle);

 // calculate E
// float Enu = (Mp*Mp) - ((Mn-Eb)*(Mn-Eb)) - (Ml*Ml) + (2.*(Mn-Eb)*El);
 float Enu = ((Mp*Mp) + (2.*(Mn-Eb)*El) - (Mn*Mn) - (Eb*Eb) + (2.*Mn*Eb)) /
              (2.*(Mn - Eb - El + (lpmom*cosbeam) ) );


 if (Enu<0.){
  float num = (Mp*Mp) - ((Mn-Eb)*(Mn-Eb)) - (Ml*Ml) + (2.*(Mn-Eb)*El);
  float denom =  (2.*(Mn - Eb - El + (lpmom*cosbeam)));
//  cout<<"num: "<<num<<endl;
//  cout<<"Ml: "<<Ml<<endl;
//  cout<<"ismu?: "<<muflg<<endl;
//  cout<<"El: "<<El<<endl;
//  cout<<"pl: "<<lpmom<<endl;
//  cout<<"denom: "<<denom<<endl;
//  cout<<"cospibm: "<<cosbeam<<endl;
 }
 //
 return Enu;
}



#endif
