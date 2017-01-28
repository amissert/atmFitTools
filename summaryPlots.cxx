#ifndef SUMMARYPLOTS_CXX
#define SUMMARYPLOTS_CXX

#include "summaryPlots.h"


/////////////////////////////////////////////////////////////
float summaryPlots::calcDipDepth(){
  
  // find the dip bin
  float pdip = 600.;
  int dipbin = pltEnuMu->FindBin(pdip);
  float dipval = pltEnuMu->GetBinContent(dipval);

  // find the max after the dip
  int maxbin = dipbin;
  float maxval = dipval;
  for (int ibin = dipbin; ibin<pltEnuMu->GetNbinsX(); ibin++){
    float binc = (float)pltEnuMu->GetBinContent(ibin);
    if (binc>maxval){
      maxval = binc;
      maxbin = ibin;
    }
  }

  // return diff
  return maxval-dipval;
}



/////////////////////////////////////////////////////////////
void summaryPlots::setPoissonErrors(TH1D* hh){

  int nbins = hh->GetNbinsX();

  for (int i=0; i<=nbins; i++){
    double content = hh->GetBinContent(i);
    hh->SetBinError(i,TMath::Sqrt(content));
  }

  return;
}



/////////////////////////////////////////////////////////////
int summaryPlots::GetCatagory(int iev, int wantnutype){
    
    // 0 -> CCQE
    // 1 -> CCnQE
    // 2 -> CCMisID
    // 3 -> NC
    // 4 -> Dead Region
   
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
  pltEnuE->Sumw2(kTRUE);
  pltEnuE->GetXaxis()->SetTitle("E_{rec} [MeV]");
  pltEnuE->SetStats(0);
  pltEnuE->SetTitle(0);
  //
  hname = nameTag.Data();
  hname.Append("_enuEBg");
  pltEnuEBg = new TH1D(hname.Data(),hname.Data(),30,0,2000);
  pltEnuEBg->Sumw2(kTRUE);
  pltEnuEBg->GetXaxis()->SetTitle("E_{rec} [MeV]");
  pltEnuEBg->SetStats(0);
  pltEnuEBg->SetTitle(0);
  //
  hname = nameTag.Data();
  hname.Append("_enuESg");
  pltEnuESg = new TH1D(hname.Data(),hname.Data(),50,0,2000);
  pltEnuESg->Sumw2(kTRUE);
  pltEnuESg->GetXaxis()->SetTitle("E_{rec} [MeV]");
  pltEnuESg->SetStats(0);
  pltEnuESg->SetTitle(0);
  //
  hname = nameTag.Data();
  hname.Append("_enuMuBg");
  pltEnuMuBg = new TH1D(hname.Data(),hname.Data(),50,0,2000);
  pltEnuMuBg->Sumw2(kTRUE);
  pltEnuMuBg->GetXaxis()->SetTitle("E_{rec} [MeV]");
  pltEnuMuBg->SetStats(0);
  pltEnuMuBg->SetTitle(0);
  //
  hname = nameTag.Data();
  hname.Append("_enuMuSg");
  pltEnuMuSg = new TH1D(hname.Data(),hname.Data(),50,0,2000);
  pltEnuMuSg->Sumw2(kTRUE);
  pltEnuMuSg->GetXaxis()->SetTitle("E_{rec} [MeV]");
  pltEnuMuSg->SetStats(0);
  pltEnuMuSg->SetTitle(0);
  //
  hname = nameTag.Data();
  hname.Append("_enuMu");
  pltEnuMu = new TH1D(hname.Data(),hname.Data(),50,0,3000);
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
    pltEnuMuCat[i] = new TH1D(hname.Data(),hname.Data(),20,0,3000);
    pltEnuMuCat[i]->GetXaxis()->SetTitle("E_{rec} [MeV]");
    pltEnuMuCat[i]->SetStats(0);
    pltEnuMuCat[i]->SetTitle(0);
    //
    hname = nameTag.Data();
    hname.Append(Form("_enuE_cat%d",i));
    pltEnuECat[i] = new TH1D(hname.Data(),hname.Data(),20,0,3000);
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
// Set up histos
void summaryPlots::clearHistos(){
 
  pltEnuE->Reset();
  pltEnuMu->Reset();
  pltEnuEBg->Reset();
  pltEnuMuBg->Reset();
  pltEnuESg->Reset();
  pltEnuMuSg->Reset();
  pltPower->Reset();
  pltSyst->Reset(); 
  pltPassMu->Reset(); 
  pltPassE->Reset(); 
  pltPower2D->Reset();
  for (int iatt=0; iatt<NATTS; iatt++){
    pltAtt[iatt]->Reset();
  }
  for (int i=0; i<NCATS; i++){
    pltEnuMuCat[i]->Reset(); 
    pltEnuECat[i]->Reset(); 
  }

  return;
  
}




/////////////////////////////////////////////////////////////////////////
// fill all histso from array
void summaryPlots::fillAllFromArray(int iev, float pow, float sys){

    // get interaction catagory
    // 0 -> CCQE
    // 1 -> CCnQE
    // 2 -> CCMisID
    // 3 -> NC
    // 4 -> Dead Region

  // get weight
  float ww = fastevents->vweight[iev];

  // for nue
  if (fastevents->vpassnue[iev]){
    int wantnutype = 12;
    // get interaction catagory
    // 0 -> CCQE
    // 1 -> CCnQE
    // 2 -> CCMisID
    // 3 -> NC
    // 4 -> Dead Region
    int icat = GetCatagory(iev, wantnutype);
    pltEnuE->Fill(fastevents->vfqenue[iev],ww);
    pltEnuECat[icat]->Fill(fastevents->vfqenue[iev],ww);
    pltPassE->Fill(fastevents->vpassnue[iev],ww);
    if (icat==0) pltEnuESg->Fill(fastevents->vfqenue[iev],ww);
    else{
      pltEnuEBg->Fill(fastevents->vfqenue[iev],ww);
    }
  }

  // for numu
  if (fastevents->vpassnumu[iev]){
    int wantnutype = 14;
    // get interaction catagory
    // 0 -> CCQE
    // 1 -> CCnQE
    // 2 -> CCMisID
    // 3 -> NC
    // 4 -> Dead Region
    int icat = GetCatagory(iev, wantnutype);
    pltEnuMu->Fill(fastevents->vfqenumu[iev],ww);
    pltEnuMuCat[icat]->Fill(fastevents->vfqenumu[iev],ww);
    pltPassMu->Fill(fastevents->vpassnumu[iev],ww);
    if (icat==0) pltEnuMuSg->Fill(fastevents->vfqenumu[iev],ww);
    else{
      pltEnuMuBg->Fill(fastevents->vfqenumu[iev],ww);
    }
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
