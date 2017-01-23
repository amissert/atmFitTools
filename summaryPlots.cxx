#ifndef SUMMARYPLOTS_CXX
#define SUMMARYPLOTS_CXX

#include "summaryPlots.h"


int summaryPlots::GetCatagory(int iev, int wantnutype){


   
    // is dead region?
    if (fastevents->vwallv[iev] < 0.) { return 4;}

    // is NC?
    if (TMath::Abs(fastevents->vmode[iev])>=30){return 3;}

    // is Mis ID?
    if (TMath::Abs(fastevents->vnutype[iev])!=wantnutype) {return 2;}

    // is CCQE?
    if (TMath::Abs(fastevents->vmode[iev])<=10){return 0;}

    // CCnQE 
    if ((TMath::Abs(fastevents->vmode[iev])>10)&&(TMath::Abs(fastevents->vmode[iev])<30)) {return 1;}
 
    cout<<"!! no catagory found"<<endl;
    
    
    return -1;

}

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
  pltEnuE = new TH1D(hname.Data(),hname.Data(),20,0,2000);
  pltEnuE->GetXaxis()->SetTitle("E_{rec} [MeV]");
  pltEnuE->SetStats(0);
  pltEnuE->SetTitle(0);
  //
  hname = nameTag.Data();
  hname.Append("_enuMu");
  pltEnuMu = new TH1D(hname.Data(),hname.Data(),20,0,5000);
  pltEnuMu->GetXaxis()->SetTitle("E_{rec} [MeV]");
  pltEnuMu->SetStats(0);
  pltEnuMu->SetTitle(0);
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

  for (int i=0; i<NCATS; i++){
    hname = nameTag.Data();
    hname.Append(Form("_enuMu_cat%d",i));
    pltEnuMuCat[i] = new TH1D(hname.Data(),hname.Data(),20,0,5000);
    pltEnuMuCat[i]->GetXaxis()->SetTitle("E_{rec} [MeV]");
    pltEnuMuCat[i]->SetStats(0);
    pltEnuMuCat[i]->SetTitle(0);
    //
    hname = nameTag.Data();
    hname.Append(Form("_enuE_cat%d",i));
    pltEnuECat[i] = new TH1D(hname.Data(),hname.Data(),20,0,1000);
    pltEnuECat[i]->GetXaxis()->SetTitle("E_{rec} [MeV]");
    pltEnuECat[i]->SetStats(0);
    pltEnuECat[i]->SetTitle(0);
  }
  
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

  // get interaction catagory
  int icat = GetCatagory(iev, fastevents->vnutype[iev]); 

  // get weight
  float ww = fastevents->vweight[iev];

  // for nue
  if (fastevents->vpassnue[iev]){
    pltEnuE->Fill(fastevents->vfqenue[iev],ww);
    pltEnuECat[icat]->Fill(fastevents->vfqenue[iev],ww);
    pltPassE->Fill(fastevents->vpassnue[iev],ww);
  }
  // for numu
  if (fastevents->vpassnumu[iev]){
    pltEnuMu->Fill(fastevents->vfqenumu[iev],ww);
    pltEnuMuCat[icat]->Fill(fastevents->vfqenumu[iev],ww);
    pltPassMu->Fill(fastevents->vpassnumu[iev],ww);
  }

  for (int i=0; i<NATTS; i++){
    pltAtt[i]->Fill(fastevents->vattribute[iev][i],ww);
  }
  if (pow>=0){
    pltPower->Fill(pow);
  }
  else{
    pltPower->Fill(fastevents->voscpower[iev][0],ww);
  }
  if (sys>=0){
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
