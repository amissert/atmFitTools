#ifndef SUMMARYPLOTS_CXX
#define SUMMARYPLOTS_CXX

#include "summaryPlots.h"

/////////////////////////////////////////////////////////////////////////
// Initialize from seed
void summaryPlots::InitToys(TH1D* hseed){
  TString hname;

  // set up toys
//  for (int ih=0; ih<NTOYPOINTS; ih++){
//    hname = nameTag.Data();
//    hname.Append(Form("_spectrum_toy%d",ih));
//    pltToySpectrum[ih] = (TH1D*)pltEnuMu->Clone(hname.Data());
    //
//    hname = nameTag.Data();
//    hname.Append(Form("_power_toy%d",ih));
//    pltToyPower[ih] = (TH1D*)pltPower->Clone(hname.Data());
    //
//    hname = nameTag.Data();
//    hname.Append(Form("_syst_toy%d",ih));
//    pltToySyst[ih] = (TH1D*)pltSyst->Clone(hname.Data());
//  }

  // set up toys
  for (int ih=0; ih<NTOYPOINTS; ih++){
    hname = nameTag.Data();
    hname.Append(Form("_spectrum_toy%d",ih));
    pltToySpectrum[ih] = (TH1D*)hseed->Clone(hname.Data());
    //
    hname = nameTag.Data();
    hname.Append(Form("_power_toy%d",ih));
    pltToyPower[ih] = (TH1D*)hseed->Clone(hname.Data());
    //
    hname = nameTag.Data();
    hname.Append(Form("_syst_toy%d",ih));
    pltToySyst[ih] = (TH1D*)hseed->Clone(hname.Data());
  }

  return;
}


/////////////////////////////////////////////////////////////////////////
// Set up histos
void summaryPlots::Init(){
  
  //
  TString hname = nameTag.Data();
  hname.Append("_enuE");
  pltEnuE = new TH1D(hname.Data(),hname.Data(),10,0,2000);
  //
  hname = nameTag.Data();
  hname.Append("_enuMu");
  pltEnuMu = new TH1D(hname.Data(),hname.Data(),10,0,2000);
  //
  hname = nameTag.Data();
  hname.Append("_Power");
  pltPower= new TH1D(hname.Data(),hname.Data(),10,0,2000);
  //
  hname = nameTag.Data();
  hname.Append("_Syst");
  pltSyst = new TH1D(hname.Data(),hname.Data(),10,0,2000);
  //
  hname = nameTag.Data();
  hname.Append("_PassNu");
  pltPassMu = new TH1D(hname.Data(),hname.Data(),100,0,3000);
  //
  hname = nameTag.Data();
  hname.Append("_PassE");
  pltPassE = new TH1D(hname.Data(),hname.Data(),100,0,3000);
  //
  for (int iatt=0; iatt<NATTS; iatt++){
    hname = nameTag.Data();
    hname.Append(Form("_att%d",iatt));
    pltAtt[iatt] = new TH1D(hname.Data(),hname.Data(),5,-1,3);
  }
  //
  hname = nameTag.Data();
  hname.Append("_pow2D");
  pltPower2D = new TH2D(hname.Data(),hname.Data(),50,0,1000,50,0,1);
 
  // set up toys
  /*
  for (int ih=0; ih<NTOYPOINTS; ih++){
    hname = nameTag.Data();
    hname.Append(Form("_spectrum_toy%d",ih));
    pltToySpectrum[ih] = (TH1D*)pltEnuMu->Clone(hname.Data());
    //
    hname = nameTag.Data();
    hname.Append(Form("_power_toy%d",ih));
    pltToyPower[ih] = (TH1D*)pltPower->Clone(hname.Data());
    //
    hname = nameTag.Data();
    hname.Append(Form("_syst_toy%d",ih));
    pltToySyst[ih] = (TH1D*)pltSyst->Clone(hname.Data());
  }
  */

  return;
  
}





/////////////////////////////////////////////////////////////////////////
// fill all histso from array
void summaryPlots::fillAllFromArray(int iev, float pow, float sys){

  //
  float ww = fastevents->vweight[iev];
  pltEnuE->Fill(fastevents->vfqenue[iev],ww);
  pltEnuMu->Fill(fastevents->vfqenue[iev],ww);
  pltPassMu->Fill(fastevents->vpassnumu[iev],ww);
  pltPassE->Fill(fastevents->vpassnue[iev],ww);
  for (int i=0; i<NATTS; i++){
    pltAtt[i]->Fill(fastevents->vattribute[iev][i],ww);
  }
  if (pow>=0){
//    plotPower2D->Fill(fastevents->vpmomv[iev],pow,ww);
    pltPower->Fill(pow);
  }
  else{
//    plotPower2D->Fill(fastevents->vpmomv[iev],fastevents->voscpower[iev],ww);
    pltPower->Fill(fastevents->voscpower[iev][0],ww);
  }
  if (sys>=0){
//    pltPower->Fill(sys);
  }
 
  return;
}

/////////////////////////////////////////////////////////////////////////
// Set pointer to large array
void summaryPlots::setLargeArray(mcLargeArray* lgarr){
  fastevents = lgarr;
}


/////////////////////////////////////////////////////////////////////////
// Constructor
summaryPlots::summaryPlots(const char* name){
  nameTag = name;
  Init();
}


#endif
