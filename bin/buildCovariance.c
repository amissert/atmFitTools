{
  gROOT->ProcessLine(".x ~/rootlogon.C");

  const int npars = 152 - 3;
  TFile *postF = new TFile("/home/xiaoyue/atmFitTools/bin/mcmctree.root", "read");
  TTree *postT = (TTree*)postF->Get("MCMCpath");
  int step;
  double par[500];
  postT->SetBranchAddress("step", &step);
  postT->SetBranchAddress("par", par);
  TFile *outF = new TFile("/home/xiaoyue/atmFitTools/bin/T2KSKcov.root", "recreate");
  TMatrixDSym *cov = new TMatrixDSym(npars);
  TMatrixDSym *cor = new TMatrixDSym(npars);
  for (int i = 0; i < npars; ++i) {
    for (int j = 0; j < npars; ++j) {
      (*cov)(i,j) = 0;
    }
  }

  int nentries = postT->GetEntries();
  double *avg = new double[npars+3];
  for (int i = 0; i < npars; ++i) { avg[i] = 0; }
  for (int i = 10000; i < nentries; ++i) {
    postT->GetEntry(i);
    for (int j = 0; j < npars+3; ++j) avg[j] += par[j];
  }
  for (int i = 0; i < npars+3; ++i) { avg[i] /= (double)(nentries - 10000); }

  for (int nn = 10000; nn < nentries; ++nn) {
    postT->GetEntry(nn);
    for (int i = 0; i < npars+3; ++i) {
      if (i > 83 && i < 87) continue;
      for (int j = 0; j < npars+3; ++j) {
	if (j > 83 && j < 87) continue;
	int ii = i; int jj = j;
	if (i > 86) ii = i-3;
	if (j > 86) jj = j-3;
	//std::cout<<i<<" "<<ii<<" "<<j<<" "<<jj<<std::endl;
	(*cov)(ii,jj) += (par[i] - avg[i]) * (par[j] - avg[j]);
      }
    }
  }
  for (int i = 0; i < npars; ++i) {
    for (int j = 0; j < npars; ++j) {
      (*cov)(i,j) /= (double)(nentries - 10000);
    }
  }
  for (int i = 0; i < npars; ++i) {
    for (int j = 0; j < npars; ++j) {
      (*cor)(i,j) = (*cov)(i,j)/TMath::Sqrt((*cov)(i,i) * (*cov)(j,j));
    }
  }

  outF->cd();
  cov->Write("covariance");
  cor->Write("correlation");
  outF->Close();
  postF->Close();
}
