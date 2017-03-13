



  // for total # of events
  for (int fvbin=0; fvbin<nfvbins; fvbin++){
    TString fname = Form("fgauss_%d",fvbin);
    float fractional_uncertainty = (hNevents[fvbin]->GetRMS()/hNevents[fvbin]->GetMean());
    float fractional_shift = TMath::Abs( hNevents[fvbin]->GetMean()- hFV[0]->GetBinContent(fvbin+1));
    fractional_shift /= hFV[0]->GetBinContent(fvbin+1);
    FVUncMap->SetBinContent(fvbin+1,fractional_uncertainty);
    FVShiftMap->SetBinContent(fvbin+1,fractional_shift);
  }

  // for CCQE
  for (int fvbin=0; fvbin<nfvbins; fvbin++){
    float fractional_uncertainty = hNeventsCCQE[fvbin]->GetRMS()/hNeventsCCQE[fvbin]->GetMean();
    float fractional_shift = TMath::Abs(hNeventsCCQE[fvbin]->GetMean() - hFVCCQE[0]->GetBinContent(fvbin+1));
    fractional_shift /= hFVCCQE[0]->GetBinContent(fvbin+1);
    FVUncMapCCQE->SetBinContent(fvbin+1,fractional_uncertainty);
    FVShiftMapCCQE->SetBinContent(fvbin+1,fractional_shift);
  }
 
  // for CCnQE
  for (int fvbin=0; fvbin<nfvbins; fvbin++){
    float fractional_uncertainty = hNeventsCCnQE[fvbin]->GetRMS()/hNeventsCCnQE[fvbin]->GetMean();
    float fractional_shift = TMath::Abs(hNeventsCCnQE[fvbin]->GetMean() - hFVCCnQE[0]->GetBinContent(fvbin+1));
    fractional_shift /= hFVCCnQE[0]->GetBinContent(fvbin+1);
    FVUncMapCCnQE->SetBinContent(fvbin+1,fractional_uncertainty);
    FVShiftMapCCnQE->SetBinContent(fvbin+1,fractional_shift);
  }

  // for CCWrong
  for (int fvbin=0; fvbin<nfvbins; fvbin++){
    float fractional_uncertainty = hNeventsCCWrong[fvbin]->GetRMS()/hNeventsCCWrong[fvbin]->GetMean();
    float fractional_shift = TMath::Abs(hNeventsCCWrong[fvbin]->GetMean() - hFVCCWrong[0]->GetBinContent(fvbin+1));
    fractional_shift /= hFVCCWrong[0]->GetBinContent(fvbin+1);
    FVUncMapCCWrong->SetBinContent(fvbin+1,fractional_uncertainty);
    FVShiftMapCCWrong->SetBinContent(fvbin+1,fractional_shift);
  }

  // for NC
  for (int fvbin=0; fvbin<nfvbins; fvbin++){
    float fractional_uncertainty = hNeventsNC[fvbin]->GetRMS()/hNeventsNC[fvbin]->GetMean();
    float fractional_shift = TMath::Abs(hNeventsNC[fvbin]->GetMean() - hFVNC[0]->GetBinContent(fvbin+1));
    fractional_shift /= hFVNC[0]->GetBinContent(fvbin+1);
    FVUncMapNC->SetBinContent(fvbin+1,fractional_uncertainty);
    FVShiftMapNC->SetBinContent(fvbin+1,fractional_shift);
  }


